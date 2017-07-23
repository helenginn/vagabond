//
//  Bucket.cpp
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Bucket.h"
#include <iostream>
#include "fftw3d.h"

/* The "checking shifts" are the shifts which need to be applied to a given
 * set of i,j,ks which are NOT within the dangerous zone of an edge - which
 * must be treated separately for speed issues */
void Bucket::makeCheckingShifts(FFTPtr map)
{
	double waterRadius = WATER_RADIUS;
	mat3x3 rebase = map->getBasisInverse();

	for (double k = -waterRadius; k < waterRadius; k += PROTEIN_SAMPLING)
	{
		for (double j = -waterRadius; j < waterRadius; j += PROTEIN_SAMPLING)
		{
			for (double i = -waterRadius; i < waterRadius; i += PROTEIN_SAMPLING)
			{
				vec3 pos = make_vec3(i, j, k);
				mat3x3_mult_vec(rebase, &pos);
				fullCheckVecs.push_back(pos);
			}
		}
	}
}

void Bucket::findBulkSolvent(FFTPtr map)
{
	makeCheckingShifts(map);
	int count = 0;
	
	for (int k = 0; k < map->nz; k++)
	{
		for (int j = 0; j < map->ny; j++)
		{
			for (int i = 0; i < map->nz; i++)
			{
				long element = map->element(i, j, k);
				MaskType type = map->getMask(element);

				if (type != MaskUnchecked)
				{
					continue;
				}

				bool foundProtein = false;

				for (int l = 0; l < fullCheckVecs.size(); l++)
				{
					int _i = fullCheckVecs[l].x + i;
					int _j = fullCheckVecs[l].y + j;
					int _k = fullCheckVecs[l].z + k;

					long adjusted = map->element(_i, _j, _k);
					

					if (map->getMask(adjusted) == MaskProtein)
					{
						foundProtein = true;
						break;
					}
				}

				if (!foundProtein)
				{
					map->setMask(element, MaskEmpty);
					count++;
				}
			}
		}
	}
}