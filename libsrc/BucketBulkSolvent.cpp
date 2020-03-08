//
//  BucketUniform.cpp
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "BucketBulkSolvent.h"
#include "Crystal.h"
#include "mat3x3.h"
#include "Molecule.h"
#include "Atom.h"
#include "Model.h"
#include <iomanip>

void BucketBulkSolvent::reportSolventContent()
{
	CrystalPtr crystal = getCrystal();
	if (crystal->isSilent())
	{
		return;
	}

	double frac = 1 - _solvent->averageAll();
	frac = (1. - frac) * 100;
	
	std::cout << std::setprecision(4) <<  "Fraction of solvent: " << frac << std::endl;
}

void BucketBulkSolvent::addSolventForConformer(int conf, int num)
{
	CrystalPtr crystal = getCrystal();

	VagFFTPtr f = crystal->getFFT();

	if (conf < 0)
	{
		_solvent = VagFFTPtr(new VagFFT(*f));
		_solvent->setAllReal(1);
		_solvent->setStatus(FFTRealSpace);
		_atomPtrs = std::vector<Atom *>(_solvent->nn(), NULL);

		for (size_t i = 0; i < crystal->atomCount(); i++)
		{
			crystal->atom(i)->addToSolventMask(_solvent, -1, &_atomPtrs, conf);
		}

		return;
	}
	
	for (size_t i = 0; i < crystal->atomCount(); i++)
	{
		AtomPtr a = crystal->atom(i);
		a->addManyToMask(_solvent, conf, num);
	}
}

void BucketBulkSolvent::adjustForVoxelVolume()
{
	mat3x3 real = _solvent->getRealBasis();
	double volume = mat3x3_volume(real);

	double conc = 0.3;
	double stuff = conc / volume;
	
	_solvent->multiplyAll(stuff);
}

void BucketBulkSolvent::addSolvent()
{
	double shrink = 0.4;
	addSolventForConformer(-1);
	_solvent->shrink(shrink);
	removeSlivers(2.0);
	setPartialStructure(_solvent);
	reportSolventContent();
	adjustForVoxelVolume();
}

void BucketBulkSolvent::removeSlivers(double maxDist)
{
	/** Maximum distance before sliver is allowed to stay in Ang */
	mat3x3 basis = _solvent->getRealBasis();
	vec3 uc_dims = empty_vec3();

	uc_dims.x = mat3x3_length(basis, 0);
	uc_dims.y = mat3x3_length(basis, 1);
	uc_dims.z = mat3x3_length(basis, 2);
	
	vec3_mult(&uc_dims, 2);

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
				_solvent->setComponent(index, 0, 0);
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
	
	for (long k = 0; k < _solvent->nz() + limits.z; k++)
	{
		for (long j = 0; j < _solvent->ny() + limits.y; j++)
		{
			for (long i = 0; i < _solvent->nx() + limits.x; i++)
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

