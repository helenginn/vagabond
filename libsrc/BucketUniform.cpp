//
//  BucketUniform.cpp
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "fftw3d.h"
#include "BucketUniform.h"

void BucketUniform::addSolvent(FFTPtr map)
{
	Bucket::findBulkSolvent(map);

	double trialDensity = 40.0;

	/* now the map's mask is set to highlight solvent */

	for (long i = 0; i < map->nn; i++)
	{
		MaskType mask = map->getMask(i);
		if (mask == MaskEmpty)
		{
			map->setElement(i, trialDensity, 0);
		}
	}
}