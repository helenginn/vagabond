//
//  BucketUniform.cpp
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "fftw3d.h"
#include "BucketBulkSolvent.h"
#include "Crystal.h"
#include "mat3x3.h"
#include "Molecule.h"

void BucketBulkSolvent::addSolvent()
{
	CrystalPtr crystal = getCrystal();

	_solvent = FFTPtr(new FFT(*crystal->getFFT()));
	mat3x3 real2frac = crystal->getReal2Frac();

	for (int i = 0; i < crystal->moleculeCount(); i++)
	{
		crystal->molecule(i)->addToMap(_solvent, real2frac, true);
	}

	_solvent->cap(1);
	_solvent->valueMinus(1);
}
