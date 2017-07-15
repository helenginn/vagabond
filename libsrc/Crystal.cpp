//
//  Crystal.cpp
//  vagabond
//
//  Created by Helen Ginn on 13/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Molecule.h"
#include "Crystal.h"
#include "fftw3d.h"
#include "vec3.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>

void Crystal::setReal2HKL(mat3x3 mat)
{
	_real2hkl = mat;
}

void Crystal::setHKL2Real(mat3x3 mat)
{
	_hkl2real = mat;
}

void Crystal::calculateMillers(FFTPtr fft)
{
	double sampling = 3.0;
	double unitScale = 1 / sampling;
	fft->setSampling(sampling);

	vec3 bounds = empty_vec3();
	bounds.x = mat3x3_length(_hkl2real, 0) * sampling;
	bounds.y = mat3x3_length(_hkl2real, 1) * sampling;
	bounds.z = mat3x3_length(_hkl2real, 2) * sampling;

	double largest = std::max(bounds.x, bounds.y);
	largest = std::max(largest, bounds.z);
	bounds.x = largest; bounds.y = largest; bounds.z = largest;
	largest *= sampling;
	largest = (int)largest;

	fft->create(bounds.x, bounds.y, bounds.z);
	fft->setScale(0, unitScale);
	fft->setScale(1, unitScale);
	fft->setScale(2, unitScale);

	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->addToMap(fft, _real2hkl);
	}

	fft->createFFTWplan(8);
	fft->fft(1);
}

void Crystal::writeCalcMillersToFile(FFTPtr fft, double resolution)
{
	double dStar = 1 / resolution;

	double aLength = mat3x3_length(_hkl2real, 0);
	double bLength = mat3x3_length(_hkl2real, 1);
	double cLength = mat3x3_length(_hkl2real, 2);

	double aLimit = aLength * dStar;
	double bLimit = bLength * dStar;
	double cLimit = cLength * dStar;

	std::ofstream file;
	file.open("tmp_name.pha");

	for (int i = -aLimit; i < aLimit; i++)
	{
		for (int j = -bLimit; j < bLimit; j++)
		{
			for (int k = -cLimit; k < cLimit; k++)
			{
				file << std::fixed << std::setprecision(1)
				<< std::setw(4) << i
				<< std::setw(4) << j
				<< std::setw(4) << k
				<< std::setw(8) << std::right << sqrt(fft->getIntensity(i, j, k))
				<<  " 1.0000  " <<
				std::setw(5) << std::right << fft->getPhase(i, j, k)
				<< std::setw(8) << 1000 << std::endl;
			}
		}
	}

	file.close();
	
}