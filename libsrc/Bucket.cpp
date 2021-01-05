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
#include <hcsrc/FileReader.h>
#include "Atom.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
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
	else
	{
		shout_at_user("Invalid choice of solvent model");
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
	
	for (int i = 0; i < 3; i++)
	{
		_averages[i] = 0;
	}
}

void Bucket::fourierTransform(int dir)
{
	/* Only care about reciprocal space */
	if (dir == 1)
	{
		_solvent->fft(FFTRealToReciprocal);
		_solvent->applySymmetry(true);

		return;
		/* F000 = nothing */
		_solvent->setComponent(0, 0, 0);
		_solvent->setComponent(0, 1, 0);
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
	return (_solvent->getReal(index) > 0.8);
}

void Bucket::abandonCalculations()
{
	return;
}

void Bucket::writeMillersToFile(std::string prefix, double maxRes)
{
	std::string solventFileOnly = prefix + "_solvent_vbond.mtz";
	CrystalPtr crystal = getCrystal();
	CSym::CCP4SPG *spg = getCrystal()->getSpaceGroup();
	
	_solvent->writeToFile(solventFileOnly, maxRes);
	
}

