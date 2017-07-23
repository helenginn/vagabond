//
//  Absolute.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Absolute.h"
#include "shared_ptrs.h"
#include "Atom.h"
#include <math.h>
#include "fftw3d.h"
#include "Element.h"
#include <iostream>
#include "Monomer.h"
#include "Polymer.h"

Absolute::Absolute(vec3 pos, double bFac, std::string element, double occValue)
{
	position = pos;
	bFactor = bFac;
	_element = element;
	_occupancy = occValue;
	_hetatm = false;
}

void Absolute::makeAtom()
{
	AtomPtr myAtom = AtomPtr(new Atom());
	myAtom->setModel(shared_from_this());
	ElementPtr element = Element::getElement(_element);
	myAtom->setElement(element);
	myAtom->setAtomName(_atomName);

	myAtom->setPosition(position);
	_atom = myAtom;
}

void Absolute::addToMolecule(MoleculePtr molecule)
{
	makeAtom();
	molecule->addAtom(_atom);

	Model::addToMolecule(molecule);
}

/*  Absolute distribution only needs to be the blurring due to the atomic
 *  B factor. The position should be provided by a different function.
 */
FFTPtr Absolute::getDistribution()
{
	FFTPtr fft = FFTPtr(new cFFTW3d());

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

void Absolute::addToMonomer(MonomerPtr monomer)
{
	makeAtom();
	monomer->addAtom(_atom);
	monomer->getPolymer()->addAtom(_atom);

	Model::addToMonomer(monomer);
}