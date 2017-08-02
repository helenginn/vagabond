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

void Crystal::realSpaceClutter()
{
	if (!_fft)
	{
		_fft = FFTPtr(new FFT());

		vec3 uc_dims = empty_vec3();
		vec3 fft_dims = empty_vec3();
		uc_dims.x = mat3x3_length(_hkl2real, 0) / PROTEIN_SAMPLING;
		uc_dims.y = mat3x3_length(_hkl2real, 1) / PROTEIN_SAMPLING;
		uc_dims.z = mat3x3_length(_hkl2real, 2) / PROTEIN_SAMPLING;

		double largest = std::max(uc_dims.x, uc_dims.y);
		largest = std::max(largest, uc_dims.z);

		fft_dims.x = largest; fft_dims.y = largest; fft_dims.z = largest;

		_fft->create(fft_dims.x, fft_dims.y, fft_dims.z);
		_fft->setupMask();

		double scaling = 1 / largest;

		_fft->setBasis(_hkl2real, scaling);
	}
	else
	{
		_fft->setAll(0);
	}


	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->addToMap(_fft, _real2frac);
	}

//	BucketPtr bucket = BucketPtr(new BucketUniform());
//	bucket->addSolvent(fft);

	_fft->createFFTWplan(8, false);
}

void Crystal::writeCalcMillersToFile(DiffractionPtr data, double resolution)
{
	if (!_fft)
	{
		shout_at_user("There is likely a bug. Cannot write\n"\
					  "calculated Miller list until it has\n"\
					  "been generated!");
	}

	realSpaceClutter();
	fourierTransform(1);
	scaleToDiffraction(data);

	double dStar = 1 / resolution;

	double aLength = mat3x3_length(_hkl2real, 0);
	double bLength = mat3x3_length(_hkl2real, 1);
	double cLength = mat3x3_length(_hkl2real, 2);

	double aLimit = aLength * dStar;
	double bLimit = bLength * dStar;
	double cLimit = cLength * dStar;

	std::string fc = _filename + "_fc.vbond.pha";
	std::ofstream fcFile;
	fcFile.open(fc);

	std::string fofc = _filename + "_fofc.vbond.pha";
	std::ofstream fofcFile;
	fofcFile.open(fofc);

	std::string twofofc = _filename + "_2fofc.vbond.pha";
	std::ofstream twofofcFile;
	twofofcFile.open(twofofc);

	/* symmetry issues */
	for (int i = -aLimit; i < aLimit; i++)
	{
		for (int j = -bLimit; j < bLimit; j++)
		{
			for (int k = 0; k < cLimit; k++)
			{
				vec3 pos = make_vec3(i, j, k);
				mat3x3_mult_vec(_real2frac, &pos);

				double intensity = _fft->getIntensity(i, j, k);
				double calcAmp = sqrt(intensity);

				double foInt = data->getFFT()->getIntensity(i, j, k);
				double foAmp = sqrt(foInt);

				if (vec3_length(pos) > dStar)
				{
					continue;
				}

				fcFile << std::fixed << std::setprecision(1)
				<< std::setw(4) << i
				<< std::setw(4) << j
				<< std::setw(4) << k
				<< std::setw(8) << std::right;

				fofcFile << std::fixed << std::setprecision(1)
				<< std::setw(4) << i
				<< std::setw(4) << j
				<< std::setw(4) << k
				<< std::setw(8) << std::right;


				twofofcFile << std::fixed << std::setprecision(1)
				<< std::setw(4) << i
				<< std::setw(4) << j
				<< std::setw(4) << k
				<< std::setw(8) << std::right;

				fcFile << calcAmp;
				fofcFile << foAmp - calcAmp;
				twofofcFile << 2 * foAmp - calcAmp;

				fcFile <<  " 1.0000  " <<
				std::setw(5) << std::right << _fft->getPhase(i, j, k)
				<< std::setw(8) << 1000 << std::endl;
				fofcFile <<  " 1.0000  " <<
				std::setw(5) << std::right << _fft->getPhase(i, j, k)
				<< std::setw(8) << 1000 << std::endl;
				twofofcFile <<  " 1.0000  " <<
				std::setw(5) << std::right << _fft->getPhase(i, j, k)
				<< std::setw(8) << 1000 << std::endl;
			}
		}
	}

	fcFile.close();
	fofcFile.close();
	twofofcFile.close();

	std::cout << "Written pha files from crystal." << std::endl;
}

double Crystal::valueWithDiffraction(DiffractionPtr data, two_dataset_op op,
									 bool verbose, double lowRes, double highRes)
{
	if (!_fft || !_fft->nn)
	{
		realSpaceClutter();
		scaleToDiffraction(data);
	}

	FFTPtr fftData = data->getFFT();
	double nLimit = std::min(fftData->nx, _fft->nx);
	nLimit /= 2;
	std::vector<double> set1, set2, free1, free2;

	double minRes = (lowRes == 0 ? 0 : 1 / lowRes);
	double maxRes = (highRes == 0 ? FLT_MAX : 1 / highRes);

	/* symmetry issues */
	for (int i = -nLimit; i < nLimit; i++)
	{
		for (int j = -nLimit; j < nLimit; j++)
		{
			for (int k = 0; k < nLimit; k++)
			{
				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(_real2frac, &ijk);
				double length = vec3_length(ijk);

				if (length < minRes || length > maxRes)
				{
					continue;
				}

				double amp1 = sqrt(fftData->getIntensity(i, j, k));
				double amp2 = sqrt(_fft->getIntensity(i, j, k));

				int isFree = fftData->getMask(i, j, k) == MaskFree;

				if (amp1 != amp1 || amp2 != amp2)
				{
					continue;
				}

				if (!isFree)
				{
					set1.push_back(amp1);
					set2.push_back(amp2);
				}
				else
				{
					free1.push_back(amp1);
					free2.push_back(amp2);
				}
			}
		}
	}

	double working = (*op)(set1, set2);

	if (verbose)
	{
		double free = (*op)(free1, free2);
		
		std::cout << "R values: " << std::setprecision(5)
		<< working;
		std::cout << ", " << free << std::endl;
	}

	return working;
}

void Crystal::applyScaleFactor(double scale, double lowRes, double highRes)
{
	double nLimit = _fft->nx;
	nLimit /= 2;
	std::vector<double> set1, set2, free1, free2;

	double minRes = (lowRes == 0 ? 0 : 1 / lowRes);
	double maxRes = (highRes == 0 ? FLT_MAX : 1 / highRes);

	/* symmetry issues */
	for (int i = -nLimit; i < nLimit; i++)
	{
		for (int j = -nLimit; j < nLimit; j++)
		{
			for (int k = 0; k < nLimit; k++)
			{
				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(_real2frac, &ijk);
				double length = vec3_length(ijk);
				long element = _fft->element(i, j, k);

				if (length < minRes || length > maxRes)
				{
					continue;
				}

				double real = _fft->getReal(element);
				double imag = _fft->getImaginary(element);

				if (real != real || imag != imag)
				{
					continue;
				}

				real *= scale;
				imag *= scale;

				_fft->setElement(element, real, imag);
			}
		}
	}
}

void Crystal::scaleToDiffraction(DiffractionPtr data)
{
	double scale = 1 / valueWithDiffraction(data, &scale_factor);
	_fft->multiplyAll(scale);

	std::vector<double> bins;
	generateResolutionBins(0, 1.0, 15, &bins);

	for (int i = 0; i < bins.size() - 1; i++)
	{
		double scale = 1 / valueWithDiffraction(data, &scale_factor, false,
										     	bins[i], bins[i + 1]);
		applyScaleFactor(scale, bins[i], bins[i + 1]);
	}
}

double Crystal::rFactorWithDiffraction(DiffractionPtr data, bool verbose)
{
	double rFactor = valueWithDiffraction(data, &r_factor, verbose);

	return rFactor;
}

void Crystal::transplantAmplitudes(DiffractionPtr data, double partsFo,
								   double partsFc)
{
	if (!_fft || !_fft->nn)
	{
		realSpaceClutter();
	}

	fourierTransform(1);

	scaleToDiffraction(data);
	rFactorWithDiffraction(data, true);

	std::cout << "Replacing Fc with Fo." << std::endl;

	FFTPtr fftData = data->getFFT();
	double nLimit = std::min(fftData->nx, _fft->nx);
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
				bool isRfree = fftData->getMask(i, j, k) == MaskFree;

				if (amp != amp || isRfree)
				{
					continue;
				}

				vec2 complex;
				long index = _fft->element(i, j, k);
				complex.x = _fft->getReal(index);
				complex.y = _fft->getImaginary(index);
				double old_amp = sqrt(complex.x * complex.x +
									  complex.y * complex.y);

				double new_amp = partsFo * amp - partsFc * old_amp;
				new_amp /= old_amp;

				complex.x *= new_amp;
				complex.y *= new_amp;

				_fft->setElement(index, complex.x, complex.y);
			}
		}
	}

	fourierTransform(-1);
}

void Crystal::tiedUpScattering()
{
	double tied = 0;
	double total = 0;

	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->tiedUpScattering(&tied, &total);
	}

	std::cout << "Tied up " << tied << " electrons out of " << total << " (";
	std::cout << 100. * sqrt(tied / total) << "%)." << std::endl;
}