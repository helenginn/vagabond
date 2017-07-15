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
#include "shared_ptrs.h"

struct vec3;

class Molecule
{
public:
	void addModel(ModelPtr model);
	void addAtom(AtomPtr atom);

	void addToMap(FFTPtr fft, mat3x3 _real2hkl);

	long int atomCount()
	{
		return atoms.size();
	}

private:
	std::vector<ModelPtr> models;
	std::vector<AtomPtr> atoms;
};

#endif /* defined(__vagabond__Molecule__) */
