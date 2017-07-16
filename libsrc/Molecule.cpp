//
//  Molecule.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "shared_ptrs.h"
#include "Molecule.h"
#include "Atom.h"
#include <float.h>
#include <iostream>
#include "fftw3d.h"
#include "mat3x3.h"

void Molecule::addModel(ModelPtr aModel)
{
	models.push_back(aModel);
}

void Molecule::addAtom(AtomPtr atom)
{
	atoms.push_back(atom);
}

void Molecule::addToMap(FFTPtr fft, mat3x3 _real2frac)
{
	for (int i = 0; i < atomCount(); i++)
	{
		atoms[i]->addToMap(fft, _real2frac);
	}
}
