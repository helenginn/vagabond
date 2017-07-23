//
//  Absolute.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Absolute.h"
#include "shared_ptrs.h"
#include "Molecule.h"
#include "Atom.h"
#include <math.h>
#include "fftw3d.h"
#include "Element.h"
#include <iostream>

Absolute::Absolute(vec3 pos, double bFac, std::string element, double occValue)
{
	position = pos;
	bFactor = bFac;
	_element = element;
	_occupancy = occValue;
}

void Absolute::addToMolecule(MoleculePtr molecule)
{
	Model::addToMolecule(molecule);

	AtomPtr myAtom = AtomPtr(new Atom());
	myAtom->addConnection(shared_from_this());
	ElementPtr element = Element::getElement(_element);
	myAtom->setElement(element);

	myAtom->setPosition(position);
	_atom = myAtom;

	molecule->addAtom(myAtom);
}

/*  Absolute distribution only needs to be the blurring due to the atomic
 *  B factor. The position should be provided by a different function.
 */
FFTPtr Absolute::getDistribution(FFTPtr *reuse)
{
	if (!(*reuse))
	{
		(*reuse) = FFTPtr(new cFFTW3d());
	}

	FFTPtr fft = *reuse;

	double n = ATOM_SAMPLING_COUNT;
	double scale = 1 / (2.0 * MAX_SCATTERING_DSTAR);
	fft->create(n);
	fft->setScales(scale);

	for (double x = -0.5; x <= 0.5; x += 1 / n)
	{
		for (double y = -0.5; y <= 0.5; y += 1 / n)
		{
			for (double z = -0.5; z <= 0.5; z += 1 / n)
			{
				double mod = MAX_SCATTERING_DSTAR * 2;
				double xAng = x * mod;
				double yAng = y * mod;
				double zAng = z * mod;

				double distSq =	(xAng * xAng + yAng * yAng + zAng * zAng);

				double value = exp((-0.25) * bFactor * distSq);
				value *= _occupancy;

				fft->setReal(x, y, z, value);
			}
		}
	}

	fft->createFFTWplan(1, false);
	
	return fft;
}