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
	_position = pos;
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

	_atom = myAtom;
}

void Absolute::addToMolecule(MoleculePtr molecule)
{
	makeAtom();
	molecule->addAtom(_atom);

	Model::addToMolecule(molecule);
}

double Absolute::getExpValue(void *object, double x, double y, double z)
{
	Absolute *me = static_cast<Absolute *>(object);
	double distSq =	(x * x + y * y + z * z);

	double bf = me->bFactor;
	double exponent = (-0.25) * bf * distSq;
	double value = exp(exponent);
	value *= me->_occupancy;

	return value;
}

/*  Absolute distribution only needs to be the blurring due to the atomic
 *  B factor. The position should be provided by a different function. */

FFTPtr Absolute::getDistribution()
{
	double n = ATOM_SAMPLING_COUNT;
	double scale = 1 / (2.0 * MAX_SCATTERING_DSTAR);

	return prepareDistribution(n, scale, this, Absolute::getExpValue);
}

void Absolute::addToMonomer(MonomerPtr monomer)
{
	makeAtom();
	monomer->addAtom(_atom);
	monomer->getPolymer()->addAtom(_atom);

	Model::addToMonomer(monomer);
}

vec3 Absolute::getPosition()
{
	return _position;
}
