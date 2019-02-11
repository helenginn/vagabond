//
//  Bucket.cpp
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Bucket.h"
#include "FileReader.h"
#include "Atom.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "fftw3d.h"
#include "RefinementNelderMead.h"
#include "Shouter.h"
#include "Crystal.h"
#include "Diffraction.h"
#include "CSV.h"
#include "Options.h"

#define CHECK_DISTANCE_STEP 0.30

Bucket::Bucket()
{
	_solvBFac = 0;
	_solvScale = 0;
	_wanted = 1;
	
	for (int i = 0; i < 3; i++)
	{
		_averages[i] = 0;
	}
}

void Bucket::scaleSolvent()
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
	fine->setEvaluationFunction(scaleSolventScore, this);
	fine->setCycles(40);
	fine->addParameter(this, getSolvScale, setSolvScale, 0.4, 0.01, "scale");
	fine->addParameter(this, getSolvBFac, setSolvBFac, 40, 1.0, "bfac");
	fine->setSilent(true);
	fine->refine();

	if (!getCrystal()->isSilent())
	{
		std::cout << "   Solvent B factor: " 
		<< getSolvBFac(this) << std::endl;
	}
	
	/** If we are doing powder analysis we don't actually want
	* 	to add the solvent */
	if (Options::shouldPowder())
	{
		return;
	}	
	
	/** Now add this into the FFT */
	
	FFTPtr fft = getCrystal()->getFFT();
	mat3x3 real2frac = getCrystal()->getReal2Frac();
	CSym::CCP4SPG *spg = getCrystal()->getSpaceGroup();
	FFTPtr fftData = _data->getFFT();

	std::vector<double> fData, fModel;
	CSVPtr csv = CSVPtr(new CSV(2, "fo", "fc"));

	vec3 nLimits = getNLimits(fftData, _solvent);

	for (int k = -nLimits.x; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.z; i < nLimits.z; i++)
			{
				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(real2frac, &ijk);
				double length = vec3_length(ijk);
				double d = 1 / length;
				double four_d_sq = (4 * d * d);
				double bFacMod = exp(-_solvBFac / four_d_sq);

				long nModel = fft->element(i, j, k);
				float realProtein = fft->data[nModel][0];
				float imagProtein = fft->data[nModel][1];
				float realSolvent = _solvent->data[nModel][0];	
				float imagSolvent = _solvent->data[nModel][1];	
				realSolvent *= _solvScale * bFacMod;
				imagSolvent *= _solvScale * bFacMod;

				float real = realProtein + realSolvent;
				float imag = imagProtein + imagSolvent;

				fft->data[nModel][0] = real;
				fft->data[nModel][1] = imag;
			}
		}
	}
}

double Bucket::scaleSolventScore(void *object)
{
	return static_cast<Bucket *>(object)->scaleAndAddSolventScore();
}

double Bucket::scaleAndAddSolventScore()
{
	FFTPtr fftData = _data->getFFT();
	FFTPtr fft = getCrystal()->getFFT();
	CSym::CCP4SPG *spg = getCrystal()->getSpaceGroup();
	double nLimit = std::min(fftData->nx, fft->nx);
	nLimit /= 2;
	mat3x3 real2frac = getCrystal()->getReal2Frac();
	mat3x3 tmp = mat3x3_transpose(real2frac);
	
	vec3 nLimits = getNLimits(fftData, _solvent);

	std::vector<double> fData, fModel;

	for (int k = -nLimits.x; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.z; i < nLimits.z; i++)
			{
				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(tmp, &ijk);

				int m, n, o;
				CSym::ccp4spg_put_in_asu(spg, i, j, k, &m, &n, &o);

				double length = vec3_length(ijk);
				double d = 1 / length;
				double four_d_sq = (4 * d * d);
				double bFacMod = exp(-_solvBFac / four_d_sq);

				int _i = 0; int _j = 0; int _k = 0;	
				CSym::ccp4spg_put_in_asu(spg, i, j, k,
				                         &_i, &_j, &_k);

				bool isRfree = (fftData->getMask(_i, _j, _k) == 0);
				if (isRfree) continue;

				float ref = sqrt(fftData->getIntensity(_i, _j, _k));

				if (ref != ref) continue;

				long nModel = fft->element(i, j, k);
				float realProtein = fft->data[nModel][0];
				float imagProtein = fft->data[nModel][1];

				float realSolvent = _solvent->data[nModel][0];	
				float imagSolvent = _solvent->data[nModel][1];	
				realSolvent *= _solvScale * bFacMod;
				imagSolvent *= _solvScale * bFacMod;

				float real = realProtein + realSolvent;
				float imag = imagProtein + imagSolvent;
				float amp = sqrt(real * real + imag * imag);

				fData.push_back(ref);
				fModel.push_back(amp);
			}
		}
	}

	double correl = correlation(fData, fModel);
	
	return -correl;
}


void Bucket::applySymOps(CSym::CCP4SPG *spaceGroup)
{
	if (spaceGroup->spg_num == 1)
	{
		return;
	}
	
	std::cout << "Solvent: ";
	_solvent->applySymmetry(spaceGroup, false);
}

void Bucket::fourierTransform(int dir)
{
	/* Only care about reciprocal space */
	if (dir == 1)
	{
		CSym::CCP4SPG *spg = getCrystal()->getSpaceGroup();

		_solvent->fft(dir);
		applySymOps(spg);
		_solvent->normalise();
		_solvent->data[0][0] = 0;
		_solvent->data[0][1] = 0;
	}
}

Atom *Bucket::nearbyAtom(int index)
{
	if (_atomPtrs.size() <= index)
	{
		return NULL;
	}

	return _atomPtrs[index];
}

bool Bucket::isSolvent(int index)
{
	return (_solvent->data[index][0] > 0.8);
}

bool Bucket::isSolvent(vec3 pos)
{
	mat3x3 real2Frac = getCrystal()->getReal2Frac();
	mat3x3_mult_vec(real2Frac, &pos);
	
	long index = _maskedRegions->elementFromFrac(pos.x, pos.y, pos.z);
	int mask = _maskedRegions->getMask(index);
	
	return (mask == 1);
}

void Bucket::processMaskedRegions()
{
	_maskedRegions = FFTPtr(new FFT(*_solvent));
	_maskedRegions->setupMask();
	FFTPtr fft = getCrystal()->getFFT();
	mat3x3 real2Frac = getCrystal()->getReal2Frac();
	mat3x3 frac2Real = mat3x3_inverse(real2Frac);
	int additions[] = {0, 0, 0};
	double sums[] = {0., 0., 0.};

	std::cout << "Total voxels: " << _maskedRegions->nn << std::endl;

	for (long i = 0; i < _maskedRegions->nn; i++)
	{
		if (_maskedRegions->data[i][0] > 0.8)
		{
			/* solvent */
			_maskedRegions->setMask(i, 1);
			additions[1]++;
			sums[1] += fft->data[i][0];
			continue;
		}
		
		additions[0]++;
		_maskedRegions->setMask(i, 0);
		sums[0] += fft->data[i][0];
		continue;

		vec3 newfrac = _maskedRegions->fracFromElement(i);
		vec3 pos = mat3x3_mult_vec(frac2Real, newfrac);
		
		AtomPtr atom = getCrystal()->getClosestAtom(pos);
		
		if (atom->isHeteroAtom())
		{
			_maskedRegions->setMask(i, 2);
			additions[2]++;
			sums[2] += fft->data[i][0];
			continue;
		}

		additions[0]++;
		sums[0] += fft->data[i][0];
	}
	
	double percentages[3] = {0., 0., 0.};
	
	for (int i = 0; i < 2; i++)
	{
		sums[i] /= (double)additions[i];
		percentages[i] = (double)additions[i] / (double)_maskedRegions->nn;
		_averages[i] = sums[i];
	}
	
	std::cout << std::setprecision(2);
	std::cout << "Protein voxels: " << additions[0];
	std::cout << " (" << percentages[0] * 100 << "%)";
	std::cout << " average value: " << std::setprecision(4)
	<< sums[0] << std::endl;

	std::cout << "Solvent voxels: " << additions[1];
	std::cout << " (" << percentages[1] * 100 << "%)";
	std::cout << " average value: " << std::setprecision(4)
	<< sums[1] << std::endl;
}

void Bucket::abandonCalculations()
{
	_solvent = FFTPtr();
}

void Bucket::writeMillersToFile(std::string prefix, double maxRes)
{
	std::string solventFileOnly = prefix + "_solvent_vbond.mtz";
	CrystalPtr crystal = getCrystal();
	std::vector<double> unitCell = crystal->getUnitCell();
	mat3x3 real2frac = getCrystal()->getReal2Frac();
	CSym::CCP4SPG *spg = getCrystal()->getSpaceGroup();
	
	_solvent->writeReciprocalToFile(solventFileOnly, maxRes,
	                                spg, unitCell, real2frac);
	
}

