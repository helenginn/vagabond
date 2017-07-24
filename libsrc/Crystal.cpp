//
//  Crystal.cpp
//  vagabond
//
//  Created by Helen Ginn on 13/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Crystal.h"
#include "fftw3d.h"
#include "vec3.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <time.h>
#include "BucketUniform.h"
#include "../libccp4/cmtzlib.h"
#include "Shouter.h"
#include "Diffraction.h"
#include "Polymer.h"

void Crystal::summary()
{
	std::cout << "|----------------" << std::endl;
	std::cout << "| Crystal summary (" << _filename << "): " << std::endl;
	std::cout << "|----------------" << std::endl;

	for (int i = 0; i < moleculeCount(); i++)
	{
		if (i > 0)
		{
			std::cout << "|-------" << std::endl;
		}
		molecule(i)->summary();
	}

	std::cout << "|----------------\n" << std::endl;
}

void Crystal::tieAtomsUp()
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->tieAtomsUp();
	}
}

void Crystal::addMolecule(MoleculePtr molecule)
{
	if (molecule->getChainID().length() <= 0)
	{
		shout_at_helen("Monomer chain ID is missing while trying\n"\
					   "to interpret PDB file.");
	}
	
	_molecules[molecule->getChainID()] = molecule;
}

void Crystal::setReal2HKL(mat3x3 mat)
{
	_real2frac = mat;
}

void Crystal::setHKL2Real(mat3x3 mat)
{
	_hkl2real = mat;
}

void Crystal::calculateMillers()
{
	fft = FFTPtr(new cFFTW3d());

	vec3 uc_dims = empty_vec3();
	vec3 fft_dims = empty_vec3();
	uc_dims.x = mat3x3_length(_hkl2real, 0) / PROTEIN_SAMPLING;
	uc_dims.y = mat3x3_length(_hkl2real, 1) / PROTEIN_SAMPLING;
	uc_dims.z = mat3x3_length(_hkl2real, 2) / PROTEIN_SAMPLING;

	double largest = std::max(uc_dims.x, uc_dims.y);
	largest = std::max(largest, uc_dims.z);

	fft_dims.x = largest; fft_dims.y = largest; fft_dims.z = largest;

	fft->create(fft_dims.x, fft_dims.y, fft_dims.z);
	fft->setupMask();

	double scaling = 1 / largest;

	fft->setBasis(_hkl2real, scaling);

	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->addToMap(fft, _real2frac);
	}

//	BucketPtr bucket = BucketPtr(new BucketUniform());
//	bucket->addSolvent(fft);

	fft->createFFTWplan(8, false);
	fft->fft(1);
	fft->multiplyAll(1e-4);
}

void Crystal::writeCalcMillersToFile(double resolution)
{
	if (!fft)
	{
		shout_at_user("There is likely a bug. Cannot write\n"\
					  "calculated Miller list until it has\n"\
					  "been generated!");
	}

	double dStar = 1 / resolution;

	double aLength = mat3x3_length(_hkl2real, 0);
	double bLength = mat3x3_length(_hkl2real, 1);
	double cLength = mat3x3_length(_hkl2real, 2);

	double aLimit = aLength * dStar;
	double bLimit = bLength * dStar;
	double cLimit = cLength * dStar;

	std::string phaFile = _filename + ".vbond.pha";
	std::ofstream file;
	file.open(phaFile);

	/* symmetry issues */
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

	std::cout << "Written " << phaFile << " from crystal." << std::endl;
}

double Crystal::valueWithDiffraction(DiffractionPtr data, two_dataset_op op,
									 bool verbose)
{
	if (!fft || !fft->nn)
	{
		calculateMillers();
	}

	FFTPtr fftData = data->getFFT();
	double nLimit = std::min(fftData->nx, fft->nx);
	nLimit /= 2;
	std::vector<double> set1, set2;

	/* symmetry issues */
	for (int i = -nLimit; i < nLimit; i++)
	{
		for (int j = -nLimit; j < nLimit; j++)
		{
			for (int k = 0; k < nLimit; k++)
			{
				double amp1 = sqrt(fftData->getIntensity(i, j, k));
				double amp2 = sqrt(fft->getIntensity(i, j, k));

				if (amp1 != amp1 || amp2 != amp2)
				{
					continue;
				}

				if (verbose)
				{
					std::cout << i << " " << j << " " << k
					<< " " << amp1 << " " << amp2 << std::endl;
				}

				set1.push_back(amp1);
				set2.push_back(amp2);
			}
		}
	}

	double value = (*op)(set1, set2);

	return value;
}

void Crystal::scaleToDiffraction(DiffractionPtr data)
{
	double scale = 1 / valueWithDiffraction(data, &scale_factor);
	fft->multiplyAll(scale);
}

double Crystal::rFactorWithDiffraction(DiffractionPtr data, bool verbose)
{
	double rFactor = valueWithDiffraction(data, &r_factor);

	if (verbose)
	{
		std::cout << "Rfactor for crystal (" << _filename << ") against data ("
		<< data->getFilename() << ") of " << std::setprecision(5)
		<< rFactor * 100 << "%." << std::endl;
	}

	return rFactor;
}

void Crystal::transplantAmplitudes(DiffractionPtr data)
{
	if (!fft || !fft->nn)
	{
		calculateMillers();
	}

	std::cout << "Replacing Fc with Fo." << std::endl;

	FFTPtr fftData = data->getFFT();
	double nLimit = std::min(fftData->nx, fft->nx);
	nLimit /= 2;
	std::vector<double> set1, set2;

	/* symmetry issues */
	for (int i = -nLimit; i < nLimit; i++)
	{
		for (int j = -nLimit; j < nLimit; j++)
		{
			for (int k = 0; k < nLimit; k++)
			{
				double amp = sqrt(fftData->getIntensity(i, j, k));

				vec2 complex;
				long index = fft->element(i, j, k);
				complex.x = fft->getReal(index);
				complex.y = fft->getImaginary(index);

				amp /= sqrt(complex.x * complex.x + complex.y * complex.y);
				complex.x *= amp;
				complex.y *= amp;

				fft->setElement(index, complex.x, complex.y);
			}
		}
	}
}