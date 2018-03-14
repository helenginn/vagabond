//
//  Distributor.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Distributor.h"
#include "fftw3d.h"
#include <iostream>

FFTPtr Distributor::prepareDistribution(double n, double scale, void *object,
                                        get_voxel_value *voxel_value)
{
	if (_precalcFFTs.count(n) && _precalcFFTs[(int)n])
	{
		_activeNum = n;
		return _precalcFFTs[n];
	}

	FFTPtr _fft = FFTPtr(new FFT());
	_fft->create(n);
	_fft->setScales(scale);
	_fft->createFFTWplan(1);
	
	for (double x = -0.5; x <= 0.5; x += 1 / n)
	{
		for (double y = -0.5; y <= 0.5; y += 1 / n)
		{
			for (double z = -0.5; z <= 0.5; z += 1 / n)
			{
				double xAng = x * _fft->getScale(0);
				double yAng = y * _fft->getScale(1);
				double zAng = z * _fft->getScale(2);

				double val = (*voxel_value)(object, xAng, yAng, zAng);

				_fft->setReal(x, y, z, val);
			}
		}
	}

	_activeNum = n;
	_precalcFFTs[n] = _fft;
	
	return _fft;
}
