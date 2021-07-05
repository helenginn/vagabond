// vagabond
// Copyright (C) 2019 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#include "FFT.h"
#include "Scaler.h"
#include "Crystal.h"
#include "Options.h"
#include "Bucket.h"
#include "Polymer.h"
#include "WeightedMap.h"

Scaler::Scaler(CrystalPtr c, DiffractionPtr d)
{
	_crystal = c;
	_data = d;
	_full = true;
	_analysis = true;
	_silent = false;
	_cycleNum = 0;
}

void Scaler::run()
{
	if (!_silent)
	{
		std::cout << "*******************************" << std::endl;
		std::cout << "\tCycle " << _cycleNum << std::endl;
	}

	prepareExplicitFFT();

	if (!_data)
	{
		_fft->writeToFile("calc_" + i_to_str(0) + ".mtz", 1.8);
		Options::flagDensityChanged();
		return;
	}

	fullSolventCalculation();

	if (_full)
	{
		writeSolvent();

		_crystal->scaleAnyPartialSet();
		scaleToData();
		tidyUp();

		makeMaps();
		_crystal->writeVagabondFile(_cycleNum);
	}
	
	if (_full && _analysis)
	{
		analyse();
	}
}

void Scaler::prepareExplicitFFT()
{
	Options::getRuntimeOptions()->disableDensityUpdate();
	_crystal->realSpaceClutter();
	_fft = _crystal->getFFT();
	_fft->fft(FFTRealToReciprocal);
	_fft->applySymmetry(false);
}

void Scaler::fullSolventCalculation()
{
	_crystal->makeBucket();
	_bucket = _crystal->getBucket();
	
	if (!_bucket)
	{
		return;
	}

	addSolvent();
	absoluteScale();
	scaleSolvent();
}

void Scaler::absoluteScale()
{
	_crystal->scaleToDiffraction(_data, false);
}

void Scaler::addSolvent()
{
	_bucket->setCrystal(_crystal);
	_bucket->addSolvent();
	_bucket->fourierTransform(1);
}

void Scaler::scaleSolvent()
{
	_bucket->setData(_data);
	_bucket->scalePartialStructure();
	_bucket->reportScale();
	_bucket->postScaleWork();
}

void Scaler::scaleToData()
{
	_crystal->scaleToDiffraction(_data, true);
}

void Scaler::writeSolvent()
{
	if (_bucket)
	{
		_bucket->writeMillersToFile("last_", 
		                            _crystal->getMaxResolution());	
	}

}

void Scaler::makeMaps()
{
	WeightedMap wMap;
	wMap.setCrystalAndData(_crystal, _data);
	wMap.createWeightedMaps();
}

void Scaler::tidyUp()
{
	_crystal->rFactorWithDiffraction(_data, true);

}

void Scaler::analyse()
{
	for (size_t i = 0; i < _crystal->moleculeCount(); i++)
	{
		MoleculePtr m = _crystal->molecule(i);
		
		if (m->isPolymer())
		{
			ToPolymerPtr(m)->ramachandranPlot();
		}
	}

}

void Scaler::findBestProbeRadius()
{
	prepareExplicitFFT();
	double bestR = FLT_MAX;
	double best_radius = 0;

	for (double r = 0; r < 0.31; r += 0.05)
	{
		Options::setProbeRadius(NULL, r);
		prepareExplicitFFT();
		fullSolventCalculation();
		scaleToData();

		double new_R = _crystal->rFactorWithDiffraction(_data, true);
		std::cout << "Probe radius " << r << " Ang; " << new_R*100
		<< "%" << std::endl;
		
		if (bestR > new_R)
		{
			bestR = new_R;
			best_radius = r;
		}
	}

	Options::setProbeRadius(NULL, best_radius);
	std::cout << std::endl;
	std::cout << "Chosen probe radius: " << best_radius << " Ang." << std::endl;
	
	run();
}

void Scaler::findProteinSampling()
{
	double bestR = FLT_MAX;
	double best_fraction = 0;

	if (Options::getProteinSampling() >= 0)
	{
		std::cout << "Protein sampling already established." << std::endl;
		return;
	}

	for (double f = 2.0; f < 4.5; f += 0.5)
	{
		Options::setSamplingFraction(NULL, f);
		_crystal->clearFFT();
		prepareExplicitFFT();
		fullSolventCalculation();
		scaleToData();

		double new_R = _crystal->rFactorWithDiffraction(_data, true);
		std::cout << "Sampling fraction 1 / " << f << "; " << new_R*100
		<< "%" << std::endl;
		
		if (bestR > new_R)
		{
			bestR = new_R;
			best_fraction = f;
		}
	}

	Options::setSamplingFraction(NULL, best_fraction);
	std::cout << std::endl;
	std::cout << "Chosen sampling fraction: 1 / " << best_fraction 
	<< "." << std::endl;
	_crystal->clearFFT();
	run();
}

void Scaler::findHetatmBSub()
{
	double bestR = FLT_MAX;
	double best_subtract = 0;

	if (Options::getBSubt() >= 0)
	{
		std::cout << "Hetatm B subtract already established." << std::endl;
		return;
	}

	for (double f = 0; f <= 20; f += 5)
	{
		Options::setBSubt(NULL, f);
		prepareExplicitFFT();
		fullSolventCalculation();
		scaleToData();

		double new_R = _crystal->rFactorWithDiffraction(_data, true);
		std::cout << "Hetatm B subtract " << f << " Ang^2; " << new_R*100
		<< "%" << std::endl;
		
		if (bestR > new_R)
		{
			bestR = new_R;
			best_subtract = f;
		}
	}

	Options::setBSubt(NULL, best_subtract);
	std::cout << std::endl;
	std::cout << "Chosen best B subtract: " << best_subtract 
	<< "." << std::endl;
	_crystal->clearFFT();
	run();
}
