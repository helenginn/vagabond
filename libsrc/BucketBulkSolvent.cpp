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
#include "Atom.h"
#include "Model.h"

void BucketBulkSolvent::addSolvent()
{
	CrystalPtr crystal = getCrystal();

	_solvent = FFTPtr(new FFT(*crystal->getFFT()));
	mat3x3 real2frac = crystal->getReal2Frac();
	_solvent->setAll(1);
	_atomPtrs = std::vector<Atom *>(_solvent->nn, NULL);
	double shrink = 0.4;

	for (size_t i = 0; i < crystal->moleculeCount(); i++)
	{
		crystal->molecule(i)->addToSolventMask(_solvent, real2frac, -1,
		                                       &_atomPtrs);
	}
	
	_solvent->shrink(shrink);
	removeSlivers();
	
	int num0 = 0;
	for (size_t i = 0; i < _solvent->nn; i++)
	{
		num0 += (_solvent->data[i][0] <= 0.5);
	}
	
	num0 *= crystal->symOpCount();

	double frac = (double)num0 / (double)_solvent->nn * 100;
	frac = 100. - frac;
	
	if (!crystal->isSilent())
	{
		std::cout << "Fraction of solvent: " << frac << std::endl;
	}
}

void BucketBulkSolvent::removeSlivers()
{
	/** Maximum distance before sliver is allowed to stay */
	double maxDist = 2.0;

	mat3x3 basis = _solvent->getBasis();
	vec3 uc_dims = empty_vec3();

	uc_dims.x = mat3x3_length(basis, 0);
	uc_dims.y = mat3x3_length(basis, 1);
	uc_dims.z = mat3x3_length(basis, 2);
	
	vec3_mult(&uc_dims, 2.0);

	uc_dims.x = (int)(maxDist / uc_dims.x + 0.5);
	uc_dims.y = (int)(maxDist / uc_dims.y + 0.5);
	uc_dims.z = (int)(maxDist / uc_dims.z + 0.5);
	
	bool changed = true;

	while (changed)
	{
		changed = sliverRemovalIteration(uc_dims);
	}
}

void BucketBulkSolvent::clearSliver(long x, long y, long z,
                                    long p, long q, long r)
{
	for (long k = z - r; k <= z + r; k++)
	{
		for (long j = y - q; j <= y + q; j++)
		{
			for (long i = x - p; i <= x + p; i++)
			{
				long index = _solvent->element(i, j, k);
				_solvent->data[index][0] = 0;
			}
		}
	}
}

bool BucketBulkSolvent::sliverRemovalIteration(vec3 limits)
{
	int lims[3];
	lims[0] = limits.x;
	lims[1] = limits.y;
	lims[2] = limits.z;
	int changed = 0;
	
	for (long k = 0; k < _solvent->nz + limits.z; k++)
	{
		for (long j = 0; j < _solvent->ny + limits.y; j++)
		{
			for (long i = 0; i < _solvent->nx + limits.x; i++)
			{
				long index = _solvent->element(i, j, k);
				float value = _solvent->getReal(index);
				
				if (value <= 0)
				{
					continue;
				}
				
				for (int n = 2; n >= 0; n--)
				{
					bool firstSide = false;

					for (int p = -lims[n]; p <= lims[n]; p++)
					{
						int px = (n == 0 ? p : 0);
						int py = (n == 1 ? p : 0);
						int pz = (n == 2 ? p : 0);
						long index = _solvent->element(i + px, j + py, k + pz);
						float value = _solvent->getReal(index);
						
						if (value <= 0 && p < 0)
						{
							firstSide = true;
						}
						else if (value <= 0 && p > 0 && firstSide)
						{
							clearSliver(i, j, k, px, py, pz);
							changed++;
							/* Make sure to trigger loop exit condition */
							p += 3 * lims[n];
							continue;
						}
						else if (p > 0 && !firstSide)
						{
							p += 3 * lims[n];
							continue;
							
							std::cout << "One side only" << std::endl;
						}
					}
				}
			}
		}
	}
	
	return changed > 0;
}

