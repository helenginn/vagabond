//
//  Crystal.h
//  vagabond
//
//  Created by Helen Ginn on 13/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Crystal__
#define __vagabond__Crystal__

#include <stdio.h>
#include <vector>
#include "shared_ptrs.h"
#include "mat3x3.h"
#include "Object.h"
#include "fftw3d.h"
#include <string>
#include "maths.h"

class Crystal : public Object
{
public:
	void addMolecule(MoleculePtr molecule)
	{
		_molecules.push_back(molecule);
	}

	long int moleculeCount()
	{
		return _molecules.size();
	}

	MoleculePtr molecule(long int i)
	{
		return _molecules[i];
	}

	void setReal2HKL(mat3x3 mat);
	void setHKL2Real(mat3x3 mat);

	void calculateMillers();
	void writeCalcMillersToFile(double resolution = 1.0);

	void scaleToDiffraction(DiffractionPtr data);
	double rFactorWithDiffraction(DiffractionPtr data);
	double valueWithDiffraction(DiffractionPtr data, two_dataset_op op,
								bool verbose = false);
	void transplantAmplitudes(DiffractionPtr data);

	void setFilename(std::string file)
	{
		_filename = file;
	}
	
private:
	std::vector<MoleculePtr> _molecules;
	std::string _filename;

	mat3x3 _hkl2real;
	mat3x3 _real2frac;

	FFTPtr fft;
};

#endif /* defined(__vagabond__Crystal__) */
