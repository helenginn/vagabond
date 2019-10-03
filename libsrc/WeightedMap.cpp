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
#include "Bucket.h"
#include "Atom.h"
#include "Options.h"
#include <iomanip>
#include "Shouter.h"
#include "Correl.h"
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
	_difft = crystal->getDifferenceMap();
}

void WeightedMap::writeCalculatedSlice()
{
	/* Back to real space */
	int nx = _fft->nx;
	FFTPtr copy = FFTPtr(new FFT(*_fft));
	copy->fft(-1);
	
	BucketPtr solv = _crystal->getBucket();
	CSVPtr calc = CSVPtr(new CSV(4, "i", "j", "d", "s"));
	
	int z = 0;

	for (int j = 0; j < _fft->ny; j++)
	{
		for (int i = 0; i < _fft->nx; i++)
		{
			int index = i + j * nx + z * nx * _fft->ny;
			Atom *atom = solv->nearbyAtom(index);
			int val = 0;
			
			if (atom)
			{
				val = 1;
			}

			if (atom && atom->getElementSymbol() == "N")
			{
				val = 2;
			}
			if (atom && atom->getElementSymbol() == "O")
			{
				val = 3;
			}
			else if (atom && atom->getElementSymbol() == "S")
			{
				val = 4;
			}
			calc->addEntry(4, (double)i, (double)j, copy->data[index][0],
			               (double)val);
		}
	}
	
	calc->writeToFile("slice_calculated.csv");
}

void WeightedMap::writeObservedSlice()
{
	CSVPtr csv = CSVPtr(new CSV(3, "i", "j", "d"));

	int nx = _fft->nx;
	long z = 0;

	for (int j = 0; j < _fft->ny; j++)
	{
		for (int i = 0; i < _fft->nx; i++)
		{
			long index = i + j * nx + z * nx * _fft->ny;
			csv->addEntry(3, (double)i, (double)j, _fft->data[index][0]);
		}
	}

	csv->writeToFile("slice_observed.csv");
	return;

	Fibonacci fib;
	fib.generateLattice(15, 3.1);
	std::vector<vec3> samples = fib.getPoints();
	mat3x3 r2f = _crystal->getReal2Frac();
	mat3x3 f2r = mat3x3_inverse(r2f);

	for (int k = 0; k < samples.size(); k++)
	{
		mat3x3_mult_vec(r2f, &samples[k]);
	}
	
	BucketPtr bucket = _crystal->getBucket();
	std::map<double, double> scores;
	
	for (double i = 0.4; i < 4.5; i += 0.4)
	{
		scores[i] = 0;
	}

	double atoms = 0;
	/*
	for (int i = 0; i < _crystal->atomCount(); i++)
	{
		AtomPtr a = _crystal->atom(i);
		if (!((a->getElementSymbol() == "O" || a->getElementSymbol() == "N")))
		{
			continue;
		}
		
		vec3 abs = a->getAbsolutePosition();
		mat3x3_mult_vec(f2r, &abs);
	*/
	
	std::cout << "Printing out solvent." << std::flush;
	double prog = 0;

	for (int k = 0; k < _fft->nz; k++)
	{
		for (int j = 0; j < _fft->ny; j++)
		{
			for (int i = 0; i < _fft->nx; i++)
			{
				long ele = _fft->element(i, j, k);
				if (bucket->nearbyAtom(ele) != NULL)
				{
					continue;
				}
				
				if (_fft->getReal(i, j, k) < 0.02)
				{
					continue;
				}
				
				vec3 frac = make_vec3((double)i / (double)_fft->nx,
				                      (double)j / (double)_fft->ny,
				                      (double)k / (double)_fft->nz);

				atoms++;

				for (std::map<double, double>::iterator it = scores.begin();
				     it != scores.end(); it++)
				{
					double dist = it->first;
					double score = FLT_MAX;
					double count = 0;

					for (int l = 0; l < samples.size(); l++)
					{
						vec3 copy = samples[l];
						vec3_mult(&copy, dist);
						vec3_add_to_vec3(&copy, frac);

						double near = _fft->getRealFromFrac(copy);
						if (near < score)
						{
							score = near;
							count++;
						}
					}

					score /= count;
					it->second += score;
				}
				
				scores[0] += _fft->getReal(i, j, k);
			}
		}
		
		prog = (double)k / _fft->nz;

		printf("\rPrinting out solvent |");
		
		for (float p = 0; p < prog; p += 0.02)
		{
			printf("=");
		}
		printf(">");

		for (float p = prog; p < 1; p += 0.02)
		{
			printf(" ");
		}
		printf("| ");

		fflush(stdout);
	}
	
	printf("\rPrinted out solvent from %.0f voxels.\n", atoms);
	
	CSVPtr solv = CSVPtr(new CSV(2, "distance", "density"));

	for (std::map<double, double>::iterator it = scores.begin();
	     it != scores.end(); it++)
	{
		it->second /= atoms;
		solv->addEntry(2, it->first, it->second);
	}
	
	solv->writeToFile("solvent_structure.csv");
}

void WeightedMap::createWeightedMaps()
{
	calculateFiguresOfMerit();
	
	Correl correl;
	correl.setCrystalAndData(_crystal, _data);
	correl.localCC();
	
	writeCalculatedSlice();

	int map = Options::getMapType();

	if (map > 0)
	{
		createVagaCoefficients();
	}
	else
	{
		create2FoFcCoefficients();
	}

	/* Back to real space */
	_crystal->fourierTransform(-1);
	_difft->fft(-1);
	
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
	
	double combined = sqrt(stdev * stdev + datadev * datadev);

	return combined;
}

double WeightedMap::oneMap(FFTPtr scratch, int slice, bool diff)
{
	FFTPtr fftData = _data->getFFT();
	double lowRes = Options::minRes();
	double minRes = (lowRes == 0 ? 0 : 1 / lowRes);
	double maxRes = _crystal->getMaxResolution(_data);
	maxRes = 1 / maxRes;
	bool ignoreRfree = Options::ignoreRFree();
	double o = -2 + (4 / MAX_SLICES) * slice;
	double weight = exp(-(o*o)/(2));

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

				long dataidx = fftData->element(_h, _k, _l);
				double fobs = fftData->data[dataidx][0];
				double sigfobs = fftData->data[dataidx][1];
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
				
				bool f000 = (i == 0 && j == 0 && k == 0);
				
				if (!f000 && ((length < minRes || length > maxRes || isAbs)
				    || (fobs != fobs || isRfree) || sigfobs != sigfobs))
				{	
					scratch->setElement(index, 0, 0);

					continue;
				}

				vec2 complex;
				complex.x = _fft->getReal(index);
				complex.y = _fft->getImaginary(index);
				double fcalc = sqrt(complex.x * complex.x +
				                      complex.y * complex.y);

				double stdev = stdevForReflection(fobs, fcalc, sigfobs,
				                                  1 / length);
				double downweight = exp(-(stdev * stdev));
				
				double phaseDev = phaseDevForWeight(downweight);
				double phi = phaseDev * o;
				
				int centric = ccp4spg_is_centric(spg, i, j, k);
				
				if (centric)
				{
					phi = 0;
				}
				
				double phase = _fft->getPhase(i, j, k);
				phi += deg2rad(phase);
				
				int shx = shellForResolution(1 / length);
				
				if (shx >= 0)
				{
					_shells[shx].count++;
					_shells[shx].phi_spread += rad2deg(phaseDev);
				}
				
				double fused = 2 * fobs - fcalc;
				
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
	FFTPtr duplicate = FFTPtr(new FFT(*_fft));
	duplicate->setAll(0);
	FFTPtr scratch = FFTPtr(new FFT(*duplicate));
	_allWeights = 0;
	
	scratch->setAll(0);
	
	for (int i = 0; i <= MAX_SLICES; i++)
	{
		double weight = oneMap(scratch, i, false);
		_allWeights += weight;
		scratch->fft(-1); /* to real space */
		scratch->multiplyAll(weight);
		FFT::addSimple(duplicate, scratch);
		scratch->setAll(0);
	}

	_difft->setAll(0);
	
	for (int i = 0; i <= MAX_SLICES; i++)
	{
		double weight = oneMap(scratch, i, true);
		scratch->fft(-1);
		scratch->multiplyAll(weight);
		FFT::addSimple(_difft, scratch);
		scratch->setAll(0);
	}
	
	double normalise = 1 / _allWeights;
	_difft->multiplyAll(normalise);
	duplicate->multiplyAll(normalise);
	
	/* To reciprocal space for writing */
	duplicate->fft(1);
	_difft->fft(1);
	
	/*
	double aveOrig = _fft->averageBoth();
	double aveDupl = duplicate->averageBoth();
	double mult = aveOrig / aveDupl;
	*/
	
	writeFile(duplicate);
	_fft->copyFrom(duplicate);
}

void WeightedMap::writeFile(FFTPtr chosen)
{
	std::string filename = _crystal->getFilename();
	double maxRes = _crystal->getMaxResolution(_data);
	CSym::CCP4SPG *spg = _crystal->getSpaceGroup();
	std::vector<double> uc = _crystal->getUnitCell();
	mat3x3 r2f = _crystal->getReal2Frac();
	
	std::string prefix = "cycle_" + i_to_str(_crystal->getCycleNum());

	std::string outputFile = prefix + "_" + filename + "_vbond.mtz";
	chosen->writeReciprocalToFile(outputFile, maxRes, spg, _data->getFFT(),
	                                 _difft, _fft);
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

				long dataidx = fftData->element(_h, _k, _l);
				double fobs = fftData->data[dataidx][0];
				double sigfobs = fftData->data[dataidx][1];
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

				if (fobs != fobs || isRfree)
				{
					_fft->setElement(index, 0, 0);
					_difft->setElement(index, 0, 0);

					continue;
				}

				vec2 complex;
				complex.x = _fft->getReal(index);
				complex.y = _fft->getImaginary(index);
				double fcalc = sqrt(complex.x * complex.x +
				                    complex.y * complex.y);

				double phase = _fft->getPhase(i, j, k);
				phase = deg2rad(phase);
				
				double stdev = stdevForReflection(fobs, fcalc, sigfobs,
				                                  1 / length);
				double downweight = exp(-(stdev * stdev));
				
				double fused = 2 * fobs - fcalc;
				double weight = exp(-stdev * stdev);
				
				if (fcalc >= fobs)
				{
//					fused = m * fobs * fobs / (d * fcalc);
				}

				complex.x = fused * cos(phase);
				complex.y = fused * sin(phase);

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
	
	writeFile(_fft);
}
