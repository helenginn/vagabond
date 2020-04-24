// force_down
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

#include "Wilson.h"
#include <libsrc/Crystal.h>
#include <libsrc/FileReader.h>
#include <libsrc/CSV.h>
#include <libsrc/maths.h>

Wilson::Wilson(VagFFTPtr vag)
{
	_fft = vag;
	_maxRes = FLT_MAX;
}

void Wilson::splitIntoShells()
{
	vec3 nLimits = getNLimits(_fft, _fft);
	int lowRes = 0;
	int ok = 0;
	double lowAng = 4.;
	CSym::CCP4SPG *spg = _fft->getSpaceGroup();

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				bool asu = CSym::ccp4spg_is_in_asu(spg, i, j, k);
				
				if (!asu)
				{
					continue;
				}

				double real = _fft->getReal(i, j, k);
				
				if (real != real)
				{
					continue;
				}

				double res = _fft->resolution(i, j, k);

				if (res < _maxRes)
				{
					_maxRes = res;
				}
				
				if (res > lowAng)
				{
					lowRes++;
				}
				
				ok++;
			}
		}
	}
	
	std::cout << "Not using " << lowRes << " reflections below "
	<< lowAng << " Å resolution." << std::endl;
	std::cout << "Maximum resolution in dataset is "
	<< _maxRes << " Å resolution." << std::endl;
	std::cout << "Using up to " << ok << " reflections for B "\
	"determination." << std::endl;
	
	int nShells = ok / 1000.;
	int resPerShell = ok / (double)nShells;
	std::cout << "Splitting into " << nShells << " shells of approx. " <<
	resPerShell << " reflections each." << std::endl;
	std::cout << std::endl;
	
	if (nShells < 6)
	{
		std::cout << " !!! ALARMINGLY LOW NUMBER OF SHELLS !!! " << std::endl;
	}

	makeShells(&_shells, lowAng, _maxRes, nShells);
}

double Wilson::getResidual(size_t maxShell)
{
	double intercept = 0;
	double gradient = 0;

	regression_line(_xs, _ys, &intercept, &gradient, maxShell);

	double total = 0;
	for (size_t i = 0; i < _xs.size() && i < maxShell; i++)
	{
		double predict = gradient * _xs[i] + intercept;
		double actual = _ys[i];
		
		double res = predict - actual;
		res *= res;
		total += res;
	}
	
	if (total != total)
	{
		return 0;
	}
	
	return total;
}

void Wilson::calculateRatios()
{
	for (size_t i = 0; i < _shells.size(); i++)
	{
		_shells[i].work1.clear();
	}

	vec3 nLimits = getNLimits(_fft, _fft);
	
	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				double res = _fft->resolution(i, j, k);
				int s = findShell(_shells, res);
				
				if (s < 0)
				{
					continue;
				}

				double real = _fft->getReal(i, j, k);
				if (real != real)
				{
					continue;
				}

				_shells[s].work1.push_back(real);
			}
		}
	}
	
	_xs.clear(); _ys.clear();
	
	for (size_t s = 0; s < _shells.size(); s++)
	{
		double ave = 0;

		for (size_t i = 0; i < _shells[s].work1.size(); i++)
		{
			ave += _shells[s].work1[i];
		}
		ave /= _shells[s].work1.size();
		double logave = log(ave);
		
		double res = (1/_shells[s].maxRes + 1/_shells[s].minRes) / 2;
		double four_dsq = 4 / (res * res);
		double x_exp = 1 / four_dsq;
		
		_xs.push_back(x_exp);
		_ys.push_back(logave);
	}
}

void Wilson::findStraightB()
{
	std::cout << "Calculating total residuals for Wilson fits for "\
	"progressively increasing numbers of shells." << std::endl;
	std::cout << std::endl;

	/*
	std::cout << "Shell (Å)\tlog(mean intensity)\tsum residual" << std::endl;
	*/
	int shellOverThreshold = -1;
	double threshold = 0.1;
	
	for (size_t i = 0; i < _xs.size(); i++)
	{
		double residual = getResidual(i);
		/*
		std::cout << pow(_xs[i], -1./3.) << "\t" << _ys[i] << "\t" <<
		residual << std::endl;
		*/
		
		if (shellOverThreshold < 0 && residual > threshold)
		{
			shellOverThreshold = i;
		}
	}
	
	if (shellOverThreshold < 0)
	{
		std::cout << "Residual never goes over " << threshold << "!" << std::endl;
		std::cout << "Nothing more to do." << std::endl;
		exit(0);
	}
	
	std::cout << "First shell over residual of " << threshold << " is at " <<
	pow(_xs[shellOverThreshold], -1./3.) << " Å." << std::endl;
	
	shellOverThreshold--;

	double intercept = 0;
	double gradient = 0;
	regression_line(_xs, _ys, &intercept, &gradient, shellOverThreshold);
	
	std::cout << "Wilson straight-B is therefore " << 
	-gradient << " Å²." << std::endl;

	_lastShell = shellOverThreshold;
	_gradient = gradient;
	_intercept = intercept;
}

void Wilson::forceDown()
{
	splitIntoShells();
	calculateRatios();
	findStraightB();
	
	std::cout << std::endl;
	std::cout << "Now correcting kinks downstream of this resolution." << 
	std::endl;

	for (size_t i = 0; i < _lastShell; i++)
	{
		_shells[i].scale = 1;
	}
	
	for (size_t i = _lastShell; i < _shells.size(); i++)
	{
		double predict = _gradient * _xs[i] + _intercept;
		double actual = _ys[i];
		
		double downweight = exp(predict) / exp(actual);
		_shells[i].scale = downweight;
	}

	std::vector<double> yTmp = _ys;
	
	correctResolutions();
	
	calculateRatios();

	CSVPtr csv = CSVPtr(new CSV(3, "res", "old", "new"));

	for (size_t i = 0; i < _xs.size(); i++)
	{
		csv->addEntry(3, _xs[i], yTmp[i], _ys[i]);
	}
	
	double yMin = FLT_MAX;
	double yMax = -FLT_MAX;
	for (size_t i = 0; i < _ys.size(); i++)
	{
		if (_ys[i] < yMin)
		{
			yMin = _ys[i];
		}

		if (_ys[i] > yMax)
		{
			yMax = _ys[i];
		}
	}
	
	std::cout << std::endl;
	
	yMin -= 0.5;
	
	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = "bfactor_plot";
	plotMap["height"] = "700";
	plotMap["width"] = "1200";
	plotMap["xHeader0"] = "res";
	plotMap["xHeader1"] = "res";
	plotMap["yHeader0"] = "new";
	plotMap["yHeader1"] = "old";
	plotMap["yMin0"] = f_to_str(yMin, 2);
	plotMap["yMin1"] = f_to_str(yMin, 2);
	plotMap["yMax0"] = f_to_str(yMax, 2);
	plotMap["yMax1"] = f_to_str(yMax, 2);

	plotMap["colour0"] = "blue";
	plotMap["colour1"] = "black";
	plotMap["xTitle0"] = "1 / (4dd)";
	plotMap["yTitle0"] = "log(val)";
	plotMap["style0"] = "line";
	plotMap["style1"] = "line";

	csv->plotPNG(plotMap);
	std::cout << "Plotting sanity check to bfactor_plot.png" << std::endl;
}

void Wilson::correctResolutions()
{
	VagFFTPtr amps = VagFFTPtr(new VagFFT(*_fft));
	amps->wipe();
	amps->setStatus(FFTReciprocalSpace);

	vec3 nLimits = getNLimits(_fft, _fft);

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				double res = _fft->resolution(i, j, k);
				int s = findShell(_shells, res);
				double real = _fft->getReal(i, j, k);
				
				if (s < 0)
				{
					amps->setReal(i, j, k, real);
					continue;
				}

				double imag = _fft->getImag(i, j, k);
				
				if (real != real)
				{
					amps->setReal(i, j, k, std::nan(""));
					continue;
				}

				real *= _shells[s].scale;
				imag *= _shells[s].scale;
				
				_fft->setReal(i, j, k, real);
				_fft->setImag(i, j, k, imag);

				amps->setReal(i, j, k, real);
			}
		}
	}
	
	std::cout << std::endl;
	std::cout << "Exporting to " << _filename << std::endl;
	std::cout << "(note only FP/SIGFP columns meaningful)" << std::endl;
	amps->writeToFile(_filename, _maxRes, _fft);
}
