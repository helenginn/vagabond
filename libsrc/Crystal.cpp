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
#include <time.h>
#include "BucketUniform.h"

void Crystal::setReal2HKL(mat3x3 mat)
{
	_real2frac = mat;
}

void Crystal::setHKL2Real(mat3x3 mat)
{
	_hkl2real = mat;
}

void Crystal::calculateMillers(FFTPtr fft)
{
	clock_t start = clock();

	double sampling = 1 / PROTEIN_SAMPLING;
	fft->setSampling(sampling);

	vec3 uc_dims = empty_vec3();
	vec3 fft_dims = empty_vec3();
	uc_dims.x = mat3x3_length(_hkl2real, 0) * sampling;
	uc_dims.y = mat3x3_length(_hkl2real, 1) * sampling;
	uc_dims.z = mat3x3_length(_hkl2real, 2) * sampling;

	double largest = std::max(uc_dims.x, uc_dims.y);
	largest = std::max(largest, uc_dims.z);

	fft_dims.x = largest; fft_dims.y = largest; fft_dims.z = largest;

	fft->create(fft_dims.x, fft_dims.y, fft_dims.z);
	fft->setupMask();

	double scaling = 1 / largest;

	fft->setMat(_hkl2real, scaling);

	clock_t end = clock();

	std::cout << "Set up map: " << (end - start) / (double)CLOCKS_PER_SEC << " seconds." << std::endl;


	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->addToMap(fft, _real2frac);
	}

	end = clock();
	std::cout << "Added atoms: " << (end - start) / (double)CLOCKS_PER_SEC << " seconds." << std::endl;

	BucketPtr bucket = BucketPtr(new BucketUniform());
	bucket->addSolvent(fft);

	end = clock();
	std::cout << "Added solvent: " << (end - start) / (double)CLOCKS_PER_SEC << " seconds." << std::endl;

//	fft->multiplyAll(0.1);
	fft->createFFTWplan(8);
	fft->fft(1);

	end = clock();

	std::cout << "Finished: " << (end - start) / (double)CLOCKS_PER_SEC << " seconds." << std::endl;

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
			for (int k = 0; k < cLimit; k++)
			{
				vec3 pos = make_vec3(i, j, k);
				mat3x3_mult_vec(_real2frac, &pos);

				if (vec3_length(pos) > dStar)
				{
					continue;
				}

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