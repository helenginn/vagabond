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

void Distributor::bFactorDistribute(FFTPtr fft, double b)
{
	mat3x3 basis = fft->getReal2Frac();
	basis = mat3x3_transpose(basis);

	for (int z = -fft->nz / 2; z < fft->nz / 2; z++)
	{
		for (int y = -fft->ny / 2; y < fft->ny / 2; y++)
		{
			for (int x = -fft->nx / 2; x < fft->nx / 2; x++)
			{
				vec3 xyz = make_vec3(x, y, z);
				mat3x3_mult_vec(basis, &xyz);
				double length = vec3_length(xyz);
				double d = 1 / length;
				double four_d_sq = (4 * d * d);
				double bFacMod = exp(-b / four_d_sq);
				int ele = fft->element(x, y, z);
				fft->data[ele][0] *= bFacMod;
				fft->data[ele][1] *= bFacMod;
			}
		}
	}
}

FFTPtr Distributor::prepareDistribution(double n, double scale, void *object,
                                        get_voxel_value *voxel_value)
{
	if (_overrideN > 0)
	{
		n = _overrideN;
	}
	
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
				double xAng = x * _fft->getScale(0) / 2;
				double yAng = y * _fft->getScale(1) / 2;
				double zAng = z * _fft->getScale(2) / 2;

				double val = (*voxel_value)(object, xAng, yAng, zAng);

				_fft->setReal(x, y, z, val);
			}
		}
	}

	_activeNum = n;
	_precalcFFTs[n] = _fft;
	
	return _fft;
}
