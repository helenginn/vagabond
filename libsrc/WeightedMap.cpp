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

#include "SpaceWarp.h"
#include "WeightedMap.h"
#include "mat3x3.h"
#include "Diffraction.h"
#include "Bucket.h"
#include "Atom.h"
#include "Options.h"
#include <iomanip>
#include "Shouter.h"
#include "CSV.h"
#include "Fibonacci.h"
#include "../libinfo/CentroidToPhase.h"

#define MAX_SLICES (25.)

WeightedMap::WeightedMap()
{

}

void WeightedMap::setCrystalAndData(CrystalPtr crystal, DiffractionPtr data)
{
	_crystal = crystal;
	_data = data;
	_fft = crystal->getFFT();
	_difft = crystal->getDiFFT();
	
	for (int i = 0; i < _fft->nn(); i++)
	{
		double real = _fft->getReal(i);
		
		if (real != real || !std::isfinite(real))
		{
			_fft->setComponent(i, 0, 0);
			_fft->setComponent(i, 1, 0);
		}
	}
}

void WeightedMap::writeCalculatedSlice()
{
	/* Back to real space */
	BucketPtr solv = _crystal->getBucket();
	VagFFTPtr mask;
	
	if (solv)
	{
		mask = solv->getSolvent();
	}
	
	if (mask)
	{
		mask->fft(FFTReciprocalToReal);
	}

	std::string cycle = "_" + i_to_str(_crystal->getCycleNum());
	mask->drawSlice(0, "solvent_slice" + cycle);

	_fft->fft(FFTReciprocalToReal);
	_fft->drawSlice(0, "calculated_slice" + cycle);
	for (int i = 0; i < _fft->nz() && Options::makeDiagnostics(); i++)
	{
//		_fft->drawSlice(i, "calculated_slice_s" + i_to_str(i));
	}

	_fft->fft(FFTRealToReciprocal);
}

void WeightedMap::writeObservedSlice()
{
	std::string cycle = "_" + i_to_str(_crystal->getCycleNum());
	
	_fft->drawSlice(0, "observed_slices" + cycle);
	for (int i = 0; i < _fft->nz() && Options::makeDiagnostics(); i++)
	{
//		_fft->drawSlice(i, "observed_slice_" + i_to_str(i));
	}
}

void WeightedMap::createWeightedMaps()
{
	calculateFiguresOfMerit();
	
	writeCalculatedSlice();

	int map = Options::getMapType();

	if (map > 0)
	{
		VagFFTPtr copy = VagFFTPtr(new VagFFT(*_fft));
		create2FoFcCoefficients(copy, true);
		copy->fft(FFTReciprocalToReal);
		double min = copy->minValue();
		double pc = -min / copy->getReal(0);
		std::cout << "Min value of F000 "; 
		std::cout << "increases electron content by " << 
		pc * 100 << "%." << std::endl;

		double f000 = _fft->getReal(0);
		f000 *= 1 + pc;
		_fft->setComponent(0, 0, f000);

		createVagaCoefficients();
	}
	else
	{
		create2FoFcCoefficients();
	}
	
	/* Back to real space */
	_crystal->fourierTransform(-1);
	_difft->fft(FFTReciprocalToReal);

	writeObservedSlice();
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
		stDiff /= aveFo;
		_shells[i].std_err = stDiff;
	}
}

int WeightedMap::shellForResolution(double res)
{
	return findShell(_shells, res);
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

	return index;
}

int centroid_qchop(double val)
{
	if (val > 1 || val < 0)
	{
		return -1;
	}
	
	int size = sizeof(centroid_to_phase_stdev) / sizeof(double) / 2;
	int min = 0;
	int max = size - 1;

	if (val <= 0) return max;
	if (val >= 1) return min;

	while (true)
	{
		int chop = (max + min) / 2;
		int higher = (val >= centroid_to_phase_stdev[chop * 2 + 1]);
		if (higher) max = chop;
		else min = chop;

		if (max - min == 1)
		{
			return max;
		}
	}
}

double WeightedMap::phaseDevForWeight(double weight)
{
	int chop = centroid_qchop(weight);

	return centroid_to_phase_stdev[chop * 2];
}

double WeightedMap::stdevForReflection(double fobs, double fcalc, 
                                       double sigfobs, double res)
{
	int shx = shellForResolution(res);
	
	double stdev = 0;

	if (shx >= 0)
	{
		stdev = _shells[shx].std_err;
	}

	stdev = 0;
	
	double datadev = sigfobs / fobs;
	double margin = (fcalc - fobs) / fobs;
	
	double combined = sqrt(stdev * stdev + 
	                       datadev * datadev +
	                       margin * margin);
	
	combined *= 2;

	return combined;
}

double WeightedMap::oneMap(VagFFTPtr scratch, int slice, bool diff)
{
	VagFFTPtr fftData = _data->getFFT();
	double lowRes = Options::minRes();
	double minRes = (lowRes == 0 ? 0 : 1 / lowRes);
	double maxRes = _crystal->getMaxResolution(_data);
	maxRes = 1 / maxRes;
	bool ignoreRfree = Options::ignoreRFree();
	double prop = -2 + (4 / MAX_SLICES) * slice;
	double weight = exp(-(prop*prop)/(2));

	if (ignoreRfree)
	{
		warn_user("You should not be ignoring Rfree.\n"\
		          "All your results and any future \n"\
		          "results derived from this are INVALID for\n"\
		          "structure determination.");
	}

	vec3 nLimits = getNLimits(_fft);
	CSym::CCP4SPG *spg = _crystal->getSpaceGroup();
	mat3x3 real2frac = _crystal->getReal2Frac();

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				long index = _fft->element(i, j, k);

				/* if no recorded observation then set element to 0 */
				if (!fftData->withinBounds(i, j, k))
				{
					scratch->setElement(index, 0, 0);
					continue;
				}
				
				int _h, _k, _l;
				CSym::ccp4spg_put_in_asu(spg, i, j, k, &_h, &_k, &_l);

				long dataidx = fftData->element(_h, _k, _l);
				double fobs = fftData->getReal(dataidx);
				double sigfobs = fftData->getImag(dataidx);

				int isAbs = CSym::ccp4spg_is_sysabs(spg, i, j, k);
				
				if (isAbs)
				{
					scratch->setElement(index, 0, 0);
					continue;
				}
				
				vec3 ijk = make_vec3(i, j, k);    
				mat3x3_mult_vec(real2frac, &ijk);
				double length = vec3_length(ijk);

				bool isFree;
				isFree = (fftData->getScratchComponent(dataidx, 0, 0) < 0.5);

				if (ignoreRfree)
				{
					isFree = 0;
				}
				
				bool f000 = (i == 0 && j == 0 && k == 0);

				if (!f000 && ((length < minRes || length > maxRes)
				    || (fobs != fobs || isFree) || sigfobs != sigfobs))
				{	
					scratch->setElement(index, 0, 0);
					continue;
				}

				vec2 complex;
				complex.x = _fft->getReal(index);
				complex.y = _fft->getImag(index);
				double fcalc = sqrt(complex.x * complex.x +
				                      complex.y * complex.y);

				double phaseDev = 0;
				if (slice == 0)
				{
					double stdev = stdevForReflection(fobs, fcalc, sigfobs,
					                                  1 / length);
					double downweight = exp(-(stdev * stdev));
					phaseDev = stdev;

					scratch->setScratchComponent(index, 0, 0, phaseDev);
				}
				else
				{
					phaseDev = scratch->getScratchComponent(index, 0, 0);
				}

				/* in radians */
				double phi = phaseDev * prop;
				
				int centric = ccp4spg_is_centric(spg, i, j, k);
				
				if (centric)
				{
					phi = 0;
				}
				
				phi = 0;
				
				double phase = _fft->getPhase(i, j, k);
				phi += deg2rad(phase);
				
				int shx = shellForResolution(1 / length);
				
				if (shx >= 0)
				{
					_shells[shx].count++;
					_shells[shx].phi_spread += rad2deg(phaseDev);
				}
				else
				{
					continue;
				}
				
				double fused = 2 * fobs - fcalc;
//				fused = fobs * exp(-(phaseDev * phaseDev));
				
				if (f000)
				{
					fused = fcalc;
				}

				if (diff)
				{
					fused = fobs - fcalc;
					if (f000) fused = 0;
				}

				complex.x = fused * cos(phi);
				complex.y = fused * sin(phi);
				
				if (complex.x != complex.x || complex.y != complex.y)
				{
					continue;
				}

				scratch->setElement(index, complex.x, complex.y);
			}
		}
	}
	
	return weight;
}

void WeightedMap::createVagaCoefficients()
{
	std::cout << "Creating Vagamap density..." << std::endl;
	/* Note: _fft currently in reciprocal space */
	VagFFTPtr duplicate = VagFFTPtr(new VagFFT(*_fft));
	duplicate->wipe();
	duplicate->setStatus(FFTRealSpace);
	VagFFTPtr scratch = VagFFTPtr(new VagFFT(*duplicate, 1));
	_allWeights = 0;
	
	scratch->makePlans();
	scratch->wipe();
	
	for (int i = 0; i <= MAX_SLICES; i++)
	{
		double weight = oneMap(scratch, i, false);
		scratch->multiplyAll(weight);
		_allWeights += weight;
		scratch->fft(FFTReciprocalToReal); /* to real space */
		duplicate->addSimple(scratch);
		scratch->multiplyAll(0);
		scratch->setStatus(FFTEmpty);
	}

	_difft->wipe();
	_difft->setStatus(FFTRealSpace);
	
	for (int i = 0; i <= MAX_SLICES; i++)
	{
		double weight = oneMap(scratch, i, true);
		scratch->fft(FFTReciprocalToReal);
		scratch->multiplyAll(weight);
		_difft->addSimple(scratch);
		scratch->wipe();
		scratch->multiplyAll(0);
		scratch->setStatus(FFTEmpty);
	}
	
	double normalise = 1 / _allWeights;
	_difft->multiplyAll(normalise);
	duplicate->multiplyAll(normalise);
	
	SpaceWarpPtr sw = Options::getActiveCrystal()->getWarp();
	sw->setCalculated(_fft);

	duplicate->fft(FFTRealToReciprocal);
	_difft->fft(FFTRealToReciprocal);

	writeFile(duplicate);
	_fft->copyFrom(duplicate);

}

void WeightedMap::writeFile(VagFFTPtr chosen)
{
	std::string filename = _crystal->getFilename();
	double maxRes = _crystal->getMaxResolution(_data);
	std::string prefix = "cycle_" + i_to_str(_crystal->getCycleNum());
	VagFFTPtr fftData = _data->getFFT();

	std::string outputFile = prefix + "_" + filename + "_vbond.mtz";
	chosen->writeToFile(outputFile, maxRes, fftData, _difft, _fft);
}

void WeightedMap::create2FoFcCoefficients(VagFFTPtr copy, bool patt)
{
//	std::cout << "Creating 2Fo-Fc density..." << std::endl;
	VagFFTPtr fftData = _data->getFFT();
	
	if (!copy)
	{
		copy = _fft;
	}

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

	vec3 nLimits = getNLimits(fftData, copy);
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

				long dataidx = fftData->element(_h, _k, _l);
				double fobs = fftData->getReal(dataidx);
				double sigfobs = fftData->getImag(dataidx);
				long index = copy->element(i, j, k);

				int isAbs = CSym::ccp4spg_is_sysabs(spg, i, j, k);
				vec3 ijk = make_vec3(i, j, k);    
				mat3x3_mult_vec(real2frac, &ijk);
				double length = vec3_length(ijk);

				bool isFree;
				isFree = (fftData->getScratchComponent(dataidx, 0, 0) < 0.5);
				
				/* If the naughty flag is set on the command line,
				 * we DO include R free in refinement - this is to show
				 * that it genuinely makes a difference */

				if (ignoreRfree)
				{
					isFree = 0;
				}
				
				if (length < minRes || length > maxRes || isAbs)
				{	
					copy->setElement(index, 0, 0);
					_difft->setElement(index, 0, 0);

					continue;
				}

				if (fobs != fobs ||
				    sigfobs != sigfobs ||
				    isFree)
				{
					copy->setElement(index, 0, 0);
					_difft->setElement(index, 0, 0);

					continue;
				}

				vec2 complex;
				complex.x = copy->getReal(index);
				complex.y = copy->getImag(index);
				double fcalc = sqrt(complex.x * complex.x +
				                    complex.y * complex.y);

				double phase = copy->getPhase(i, j, k);
				phase = deg2rad(phase);
				
				double stdev = stdevForReflection(fobs, fcalc, sigfobs,
				                                  1 / length);
				double downweight = exp(-(stdev * stdev));
				
				double fused = 2 * fobs - fcalc;
				
				complex.x = fused * cos(phase);
				complex.y = fused * sin(phase);
				
				if (patt)
				{
					complex.x = fobs;
					complex.y = 0;
				}

				if (complex.x != complex.x || complex.y != complex.y)
				{
					continue;
				}

				copy->setElement(index, complex.x, complex.y);
				
				bool f000 = (i == 0 && j == 0 && k == 0);
				
				if (f000)
				{
					continue;
				}
				
				fused = fobs - fcalc;

				complex.x = fused * cos(phase);
				complex.y = fused * sin(phase);

				if (complex.x != complex.x || complex.y != complex.y)
				{
					continue;
				}

				_difft->setElement(index, complex.x, complex.y);
			}
		}
	}
	
	writeFile(copy);
}
