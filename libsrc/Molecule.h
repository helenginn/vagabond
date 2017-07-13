//
//  Molecule.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Molecule__
#define __vagabond__Molecule__

#include <stdio.h>
#include <vector>
#include "mat3x3.h"

struct vec3;

class Molecule
{
public:
	void addModel(ModelPtr model);
	void addAtom(AtomPtr atom);

	void setReal2HKL(mat3x3 mat);
	void setHKL2Real(mat3x3 mat);

	void calculateMillers(FFTPtr fft);

	long int atomCount()
	{
		return atoms.size();
	}

private:
	mat3x3 _hkl2real;
	mat3x3 _real2hkl;
	std::vector<ModelPtr> models;
	std::vector<AtomPtr> atoms;
};

#endif /* defined(__vagabond__Molecule__) */
