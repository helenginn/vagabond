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
#include <map>

typedef double(get_voxel_value)(void *obj, double x, double y, double z);

class Distributor
{
public:
	Distributor()
	{
		_calculated = false;
		_activeNum = 0;
	}

	virtual FFTPtr getDistribution(bool quick = false, int n = -1) = 0;

	void recalculate()
	{
		_precalcFFTs.clear();
	}
protected:
	bool _calculated;

	FFTPtr getDistributionCopy()
	{
		if (_activeNum <= 0) return FFTPtr();

		FFTPtr newPtr;
		newPtr.reset(new FFT(*_precalcFFTs[_activeNum]));
		return newPtr;
	}

	FFTPtr prepareDistribution(double n, double scale, void *object,
	                           get_voxel_value *voxel_value);

	FFTPtr getFFT(int n)
	{
		return _precalcFFTs[n];	
	}
	
	virtual std::string getClassName()
	{
		return "Distributor";
	}

private:
	int _activeNum;
	std::map<int, FFTPtr> _precalcFFTs;
};

#endif /* defined(__vagabond__Distributor__) */
