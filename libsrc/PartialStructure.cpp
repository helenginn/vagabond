// Vagabond
// Copyright (C) 2017-2018 Helen Ginn
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

#include "PartialStructure.h"
#include "maths.h"
#include "fftw3d.h"
#include "Diffraction.h"
#include "CSV.h"
#include "RefinementNelderMead.h"
#include "Options.h"
#include "Crystal.h"
#include "Shouter.h"

void PartialStructure::setStructure(FFTPtr refPart)
{
	FFTPtr fft = getCrystal()->getFFT();
	CSym::CCP4SPG *spg = getCrystal()->getSpaceGroup();

	_partial = FFTPtr(new FFT(*fft));
	_partial->setAll(0);

	vec3 nLimits = getNLimits(_partial, refPart);

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				int _h, _k, _l;
				int sym = CSym::ccp4spg_put_in_asu(spg, i, j, k, 
				                                   &_h, &_k, &_l);

				long index = _partial->element(i, j, k);
				long asuIndex = _partial->element(_h, _k, _l);
				long pIndex = refPart->element(_h, _k, _l);
				
				double x = refPart->data[pIndex][0];
				double y = refPart->data[pIndex][1];
				
				/* establish amp & phase */
				double amp = sqrt(x * x + y * y);
				double myPhase = atan2(y, x);

				if (sym % 2 == 0)
				{
					myPhase *= -1;
				}
				
				int symop = (sym - 1) / 2;

				/* calculate phase shift for symmetry operator */
				float *trn = spg->symop[symop].trn;
				
				/* rotation */
				double shift = (float)i * trn[0];
				shift += (float)j * trn[1];
				shift += (float)k * trn[2];
				
				shift = fmod(shift, 1.);

				/*  apply shift when filling in other sym units */
				double newPhase = myPhase + shift * 2 * M_PI;

				/* debug */
				float symx = fft->data[index][0];
				float symy = fft->data[index][1];
				double orig = atan2(symy, symx);
				
				x = amp * cos(newPhase);
				y = amp * sin(newPhase);

				_partial->data[index][0] = x;
				_partial->data[index][1] = y;
			}
		}
	}
}

void PartialStructure::reportScale()
{
	if (!getCrystal()->isSilent())
	{
		std::cout << "   Partial scale: " << getSolvScale(this) << ", "
		"B factor: " << getSolvBFac(this) << std::endl;
	}
}

void PartialStructure::scalePartialStructure()
{
	if (!_data)
	{
		shout_at_helen("Need diffraction data to scale solvent");
	}
	
	setSolvScale(this, 0.0);
	setSolvBFac(this, 40);
	
	NelderMeadPtr fine;
	fine = NelderMeadPtr(new RefinementNelderMead());
	fine->setJobName("solvent_scale_fine");
	fine->setEvaluationFunction(scalePartialScore, this);
	fine->setCycles(40);
	fine->addParameter(this, getSolvScale, setSolvScale, 0.4, 0.01, "scale");
	fine->addParameter(this, getSolvBFac, setSolvBFac, 40, 1.0, "bfac");
	fine->setSilent(true);
	fine->refine();

	/** If we are doing powder analysis we don't actually want
	* 	to add the solvent */
	if (Options::shouldPowder())
	{
		return;
	}	
	
	/** Now add this into the FFT */
	
	FFTPtr fft = getCrystal()->getFFT();
	mat3x3 real2frac = getCrystal()->getReal2Frac();
	real2frac = mat3x3_transpose(real2frac);
	CSym::CCP4SPG *spg = getCrystal()->getSpaceGroup();
	FFTPtr fftData = _data->getFFT();

	std::vector<double> fData, fModel;
	CSVPtr csv = CSVPtr(new CSV(2, "fo", "fc"));

	vec3 nLimits = getNLimits(fftData, _partial);

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				long nModel = fft->element(i, j, k);
				float realProtein = fft->data[nModel][0];
				float imagProtein = fft->data[nModel][1];
				float realPartial = _partial->data[nModel][0];	
				float imagPartial = _partial->data[nModel][1];	
				
				if (realProtein != realProtein ||
				    realPartial != realPartial)
				{
					continue;
				}

				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(real2frac, &ijk);
				double length = vec3_length(ijk);
				double d = 1 / length;
				double four_d_sq = (4 * d * d);
				double bFacMod = exp(-_solvBFac / four_d_sq);

				realPartial *= _solvScale * bFacMod;
				imagPartial *= _solvScale * bFacMod;

				float real = realProtein + realPartial;
				float imag = imagProtein + imagPartial;
				
				if (real != real || imag != imag)
				{
					continue;
				}

				fft->data[nModel][0] = real;
				fft->data[nModel][1] = imag;
			}
		}
	}
}

double PartialStructure::scalePartialScore(void *object)
{
	return static_cast<PartialStructure *>(object)->scaleAndAddPartialScore();
}

double PartialStructure::scaleAndAddPartialScore()
{
	FFTPtr fftData = _data->getFFT();
	FFTPtr fft = getCrystal()->getFFT();
	CSym::CCP4SPG *spg = getCrystal()->getSpaceGroup();
	mat3x3 real2frac = getCrystal()->getReal2Frac();
	mat3x3 tmp = mat3x3_transpose(real2frac);
	
	vec3 nLimits = getNLimits(fftData, _partial);

	std::vector<double> fData, fModel;
	
	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				long nModel = fft->element(i, j, k);
				float realProtein = fft->data[nModel][0];
				float imagProtein = fft->data[nModel][1];
				float realPartial = _partial->data[nModel][0];	
				float imagPartial = _partial->data[nModel][1];	
				
				if (realProtein != realProtein ||
				    realPartial != realPartial)
				{
					continue;
				}

				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(tmp, &ijk);

				int m, n, o;
				int asu = CSym::ccp4spg_put_in_asu(spg, i, j, k, &m, &n, &o);
				
				double length = vec3_length(ijk);
				double d = 1 / length;
				double four_d_sq = (4 * d * d);
				double bFacMod = exp(-_solvBFac / four_d_sq);

				bool isRfree = (fftData->getMask(m, n, o) == 0);
				if (isRfree) continue;

				float ref = sqrt(fftData->getIntensity(m, n, o));

				if (ref != ref) continue;

				realPartial *= _solvScale * bFacMod;
				imagPartial *= _solvScale * bFacMod;
				
				float real = realProtein + realPartial;
				float imag = imagProtein + imagPartial;
				float amp = sqrt(real * real + imag * imag);
				
				fData.push_back(ref);
				fModel.push_back(amp);
			}
		}
	}


	double correl = correlation(fData, fModel);
	
	return -correl;
}
