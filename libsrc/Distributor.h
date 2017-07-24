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
#include "fftw3d.h"

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
	bool _calculated;

	FFTPtr getDistributionCopy()
	{
		return std::make_shared<cFFTW3d>(*_fft);
	}

	FFTPtr prepareDistribution(double n, double scale, void *object,
							   get_voxel_value *voxel_value);

private:
	FFTPtr _fft;
};

#endif /* defined(__vagabond__Distributor__) */
