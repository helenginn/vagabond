//
//  Distributor.h
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Distributor__
#define __vagabond__Distributor__

#include <stdio.h>
#include "shared_ptrs.h"

typedef double(get_voxel_value)(void *obj, double x, double y, double z);

class Distributor
{
public:
	Distributor()
	{
		_calculated = false;
	}
	
	virtual FFTPtr getDistribution() = 0;

protected:
	FFTPtr _fft;
	bool _calculated;

	FFTPtr prepareDistribution(get_voxel_value *func,
							   bool fftNow = false);
};

#endif /* defined(__vagabond__Distributor__) */
