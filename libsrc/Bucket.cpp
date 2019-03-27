//
//  Bucket.cpp
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Bucket.h"
#include "BucketBulkSolvent.h"
#include "BucketPerStrand.h"
#include "FileReader.h"
#include "Atom.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "fftw3d.h"
#include "Shouter.h"
#include "Crystal.h"
#include "Diffraction.h"
#include "CSV.h"
#include "Options.h"

#define CHECK_DISTANCE_STEP 0.30

BucketPtr Bucket::chosenBucket()
{
	BucketPtr bucket;

	if (Options::getAddSolvent() == 1)
	{
		bucket = BucketPtr(new BucketBulkSolvent());
	}
	else if (Options::getAddSolvent() == 2)
	{
		bucket = BucketPtr(new BucketPerStrand());
	}
	
	return bucket;
}

void Bucket::reportScale()
{
	if (!getCrystal()->isSilent())
	{
		std::cout << "   Solvent scale: " << getSolvScale(this) << ", "
		"B factor: " << getSolvBFac(this) << std::endl;
	}
}


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
	return;
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

