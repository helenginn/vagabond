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

#include "WeightedMap.h"
#include "mat3x3.h"
#include "Diffraction.h"
#include "Options.h"
#include "fftw3d.h"
#include "Shouter.h"

WeightedMap::WeightedMap()
{

}

void WeightedMap::setCrystalAndData(CrystalPtr crystal, DiffractionPtr data)
{
	_crystal = crystal;
	_data = data;
	_fft = crystal->getFFT();
	_difft = crystal->getDifferenceMap();
}

void WeightedMap::createWeightedMaps()
{
	createWeightedCoefficients();

	/* Back to real space */
	_crystal->fourierTransform(-1);
	_difft->fft(-1);
}

void WeightedMap::calculateFiguresOfMerit()
{
	_shells = _crystal->getShells();
	
	
	double sumFo = 0;
	double sumFc = 0;

	for (int i = 0; i < _shells.size(); i++)
	{
		for (int j = 0; j < _shells[i].work1.size(); j++)
		{
			double fo = _shells[i].work1[j];
			double fc = _shells[i].work2[j];

			sumFo += fo;
			sumFc += fc;
		}

		double aveFo = sumFo / (double)_shells[i].work1.size();
		double aveFc = sumFc / (double)_shells[i].work1.size();
		double sumDiff = 0;

		std::cout << "Average Fo: " << aveFo << std::endl;

		for (int j = 0; j < _shells[i].work1.size(); j++)
		{
			double fo = _shells[i].work1[j];
			double fc = _shells[i].work2[j];

			double diff = pow(fo - fc, 2) / (fo * aveFo);
			sumDiff += diff;
		}

		double aveDiff = sumDiff / (double)_shells[i].work1.size();
		std::cout << _shells[i].maxRes << " " << aveDiff << std::endl;
	}
}

void WeightedMap::create2FoFcCoefficients()
{
	FFTPtr fftData = _data->getFFT();

	double lowRes = Options::minRes();
	double minRes = (lowRes == 0 ? 0 : 1 / lowRes);
	double maxRes = _crystal->getMaxResolution(_data);
	maxRes = 1 / maxRes;
	bool ignoreRfree = Options::ignoreRFree();
	double partsFo = 2;
	double partsFc = 1;
	
	if (ignoreRfree)
	{
		warn_user("You should not be ignoring Rfree.\n"\
		          "All your results and any future \n"\
		          "results derived from this are INVALID for\n"\
		          "structure determination.");
	}
	
	std::cout << "Mixing in data for weighed density map..." << std::endl;

	vec3 nLimits = getNLimits(fftData, _fft);
	CSym::CCP4SPG *spg = _crystal->getSpaceGroup();
	mat3x3 real2frac = _crystal->getReal2Frac();

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				int _h, _k, _l;
				CSym::ccp4spg_put_in_asu(spg, i, j, k, &_h, &_k, &_l);

				double obs_amp = sqrt(fftData->getIntensity(_h, _k, _l));
				long index = _fft->element(i, j, k);

				int isAbs = CSym::ccp4spg_is_sysabs(spg, i, j, k);
				vec3 ijk = make_vec3(i, j, k);    
				mat3x3_mult_vec(real2frac, &ijk);
				double length = vec3_length(ijk);

				bool isRfree = (fftData->getMask(_h, _k, _l) == 0);
				
				/* If the naughty flag is set on the command line,
				 * we DO include R free in refinement - this is to show
				 * that it genuinely makes a difference */

				if (ignoreRfree)
				{
					isRfree = 0;
				}
				
				if (length < minRes || length > maxRes || isAbs)
				{	
					_fft->setElement(index, 0, 0);
					_difft->setElement(index, 0, 0);

					continue;
				}

				if (obs_amp != obs_amp || isRfree)
				{
					_fft->setElement(index, 0, 0);
					_difft->setElement(index, 0, 0);

					continue;
				}


				vec2 complex;
				complex.x = _fft->getReal(index);
				complex.y = _fft->getImaginary(index);
				double calc_amp = sqrt(complex.x * complex.x +
				                      complex.y * complex.y);

				double new_amp = calc_amp;
				new_amp = partsFo * obs_amp - partsFc * calc_amp;
				double rescale = new_amp / calc_amp;

				double diff_scale = obs_amp - calc_amp;
				diff_scale /= calc_amp;

				vec2 diff_complex = complex;
				complex.x *= rescale;
				complex.y *= rescale;

				if (complex.x != complex.x || complex.y != complex.y)
				{
					continue;
				}

				_fft->setElement(index, complex.x, complex.y);
				
				bool f000 = (i == 0 && j == 0 && k == 0);
				
				if (f000)
				{
					continue;
				}
				
				diff_complex.x *= diff_scale;
				diff_complex.y *= diff_scale;

				if (diff_complex.x != diff_complex.x || 
				    diff_complex.y != diff_complex.y)
				{
					continue;
				}

				_difft->setElement(index, diff_complex.x, diff_complex.y);
			}
		}
	}
}
