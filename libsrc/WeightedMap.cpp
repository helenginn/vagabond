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

#define MAX_SLICES (25)

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
	calculateFiguresOfMerit();
//	create2FoFcCoefficients();
	createVagaCoefficients();

	/* Back to real space */
	_crystal->fourierTransform(-1);
	_difft->fft(-1);
}

void WeightedMap::calculateFiguresOfMerit()
{
	_shells = _crystal->getShells();
	
	for (int i = 0; i < _shells.size(); i++)
	{
		double sumFo = 0;
		double sumFc = 0;

		for (int j = 0; j < _shells[i].work1.size(); j++)
		{
			double fo = _shells[i].work1[j];
			double fc = _shells[i].work2[j];

			sumFo += fo;
			sumFc += fc;
		}

		double aveFo = sumFo / (double)_shells[i].work1.size();
		double aveFc = sumFc / (double)_shells[i].work1.size();
		double scale = aveFo / aveFc;
		double sumDiff = 0;
		_shells[i].aveFo = aveFo;
		
		for (int j = 0; j < _shells[i].work1.size(); j++)
		{
			double fo = _shells[i].work1[j];
			double fc = _shells[i].work2[j];

			double diff = fabs(fo - scale * fc);
			diff *= diff;
			sumDiff += diff;
		}

		double stDiff = sqrt(sumDiff / _shells[i].work1.size());
		_shells[i].std_err = stDiff / aveFo;
	}
}

/* y between -1 and 1 */
double evaluationOnCircle(double fobs, double fcalc, double stdev, double y)
{
	double d = fobs * fobs + fcalc * fcalc;
	double n = -2 * fobs * fcalc;

	double expo = d + n * y;
	double stdexpo = d + n;
	expo *= -1 / stdev;
	stdexpo *= -1 / stdev;
	
	return exp(expo - stdexpo);
}

double cutoff_position(double fobs, double fcalc, double stdev, double cut)
{
	double d = fobs * fobs + fcalc * fcalc;
	double n = -2 * fobs * fcalc;
	d /= -stdev; n /= -stdev;
	
	double y = (log(cut) + n) / n;
	
	if (y < -1)
	{
		y = -1;
	}
	
	return y;
}

double WeightedMap::stdevForReflection(double fobs, double fcalc, double res)
{
	int index = -1;

	for (int l = 0; l < _shells.size(); l++)
	{
		if (res <= _shells[l].minRes &&
		    res > _shells[l].maxRes)
		{
			index = l;
			break;
		}
	}

	if (index < 0)
	{
		return HUGE_VAL;
	}

	double std_error = _shells[index].std_err;
	double aveFo = _shells[index].aveFo;
	double std_stdev = std_error * aveFo;
	double my_dev = fabs(fobs - fcalc);
	double combined_dev = sqrt(my_dev * my_dev + std_stdev * std_stdev);
	double stdev = combined_dev;

	return stdev;
}

void WeightedMap::oneMap(FFTPtr scratch, int slice, bool diff)
{
	FFTPtr fftData = _data->getFFT();
	double divide = (MAX_SLICES - 1) / 2;
	double lowRes = Options::minRes();
	double minRes = (lowRes == 0 ? 0 : 1 / lowRes);
	double maxRes = _crystal->getMaxResolution(_data);
	maxRes = 1 / maxRes;
	bool ignoreRfree = Options::ignoreRFree();

	if (ignoreRfree)
	{
		warn_user("You should not be ignoring Rfree.\n"\
		          "All your results and any future \n"\
		          "results derived from this are INVALID for\n"\
		          "structure determination.");
	}

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

				double fobs = sqrt(fftData->getIntensity(_h, _k, _l));
				long index = _fft->element(i, j, k);

				int isAbs = CSym::ccp4spg_is_sysabs(spg, i, j, k);
				vec3 ijk = make_vec3(i, j, k);    
				mat3x3_mult_vec(real2frac, &ijk);
				double length = vec3_length(ijk);

				bool isRfree = (fftData->getMask(_h, _k, _l) == 0);

				if (ignoreRfree)
				{
					isRfree = 0;
				}
				
				if ((length < minRes || length > maxRes || isAbs)
				    || (fobs != fobs || isRfree))
				{	
					scratch->setElement(index, 0, 0);

					continue;
				}

				vec2 complex;
				complex.x = _fft->getReal(index);
				complex.y = _fft->getImaginary(index);
				double fcalc = sqrt(complex.x * complex.x +
				                      complex.y * complex.y);

				double stdev = stdevForReflection(fobs, fcalc, 1 / length);
//				stdev *= stdev;
				double yCut = cutoff_position(fobs, fcalc, stdev, 0.05);
				double portion = (1 - yCut) / divide;
				double yCurr = yCut + portion * slice;
				if (yCurr > 1) 
				{
					yCurr = 2 - yCurr;
				}

				double phi = acos(yCurr);
				
				if (slice < divide) phi *= -1;
				
				double phase = _fft->getPhase(i, j, k);
				phi += deg2rad(phase);
				double weight = evaluationOnCircle(fobs, fcalc, stdev, yCurr);

				if (diff)
				{
					complex.x = (fobs - fcalc) * weight * cos(phi);
					complex.y = (fobs - fcalc) * weight * sin(phi);
				}
				else
				{
					complex.x = fobs * weight * cos(phi);
					complex.y = fobs * weight * sin(phi);
				}
				
				if (complex.x != complex.x || complex.y != complex.y)
				{
					continue;
				}

				scratch->setElement(index, complex.x, complex.y);
				
				bool f000 = (i == 0 && j == 0 && k == 0);
				
				if (f000)
				{
					continue;
				}
			}
		}
	}
}

void WeightedMap::createVagaCoefficients()
{
	std::cout << "Creating Vagamap density..." << std::endl;
	FFTPtr duplicate = FFTPtr(new FFT(*_fft));
	FFTPtr scratch = FFTPtr(new FFT(*_fft));
	
	for (int i = 0; i < MAX_SLICES; i++)
	{
		oneMap(scratch, i, false);
		scratch->fft(-1);
		FFT::addSimple(duplicate, scratch);
		scratch->setAll(0);
	}
	

	_difft->setAll(0);
	
	for (int i = 0; i < MAX_SLICES; i++)
	{
		oneMap(scratch, i, true);
		scratch->fft(-1);
		FFT::addSimple(_difft, scratch);
		scratch->setAll(0);
	}
	
	_difft->multiplyAll(0.5);
	FFT::addSimple(duplicate, _difft);
	_difft->multiplyAll(2.0);
	
	duplicate->fft(1);
	_fft->copyFrom(duplicate);
	_difft->fft(1);
}

void WeightedMap::create2FoFcCoefficients()
{
	std::cout << "Creating 2Fo-Fc density..." << std::endl;
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
