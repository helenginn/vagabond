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

class Crystal
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

	void calculateMillers(FFTPtr fft);
	void writeCalcMillersToFile(FFTPtr fft, double resolution = 1.0);

private:
	std::vector<MoleculePtr> _molecules;

	mat3x3 _hkl2real;
	mat3x3 _real2frac;

};

#endif /* defined(__vagabond__Crystal__) */
