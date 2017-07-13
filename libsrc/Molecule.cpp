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

void Molecule::calculateMillers(FFTPtr fft)
{
	double sampling = 3.0;
	fft->setSampling(sampling);

	vec3 bounds = empty_vec3();
	bounds.x = mat3x3_length(_hkl2real, 0);
	bounds.y = mat3x3_length(_hkl2real, 1);
	bounds.z = mat3x3_length(_hkl2real, 2);

	double largest = std::max(bounds.x, bounds.y);
	largest = std::max(largest, bounds.z);
	largest *= sampling;

	fft->create(largest, largest, largest);
	fft->setScale(0, largest);
	fft->setScale(1, largest);
	fft->setScale(2, largest);

	mat3x3 sample2hkl = _real2hkl;
	mat3x3_scale(&sample2hkl, 1 / sampling, 1 / sampling, 1 / sampling);

	for (int i = 0; i < atomCount(); i++)
	{
		atoms[i]->addToMap(fft, sample2hkl);
	}
}

void Molecule::setReal2HKL(mat3x3 mat)
{
	_real2hkl = mat;
}

void Molecule::setHKL2Real(mat3x3 mat)
{
	_hkl2real = mat;
}
