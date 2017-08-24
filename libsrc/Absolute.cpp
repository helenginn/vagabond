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
#include "maths.h"

Absolute::Absolute(vec3 pos, double bFac, std::string element, double occValue)
{
	_position = pos;
	_bFactor = bFac;
	_element = element;
	trim(_element);
	_occupancy = occValue;
	_hetatm = false;
	_usingTensor = false;
	_tensor = make_mat3x3();
}

void Absolute::makeAtom()
{
	AtomPtr myAtom = AtomPtr(new Atom());
	myAtom->setModel(shared_from_this());
	myAtom->setInitialPosition(_position);
	myAtom->setInitialBFactor(_bFactor);
	myAtom->setAtomNum(_atomNum);
	ElementPtr element = Element::getElement(_element);
	myAtom->setElement(element);
	myAtom->setAtomName(_atomName);
	myAtom->findAtomType(_resName);

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
	double aniso = 0;

	if (me->_usingTensor)
	{
		vec3 recipVec = make_vec3(x, y, z);
		mat3x3_mult_vec(me->_tensor, &recipVec);
		recipVec.x *= x;
		recipVec.y *= y;
		recipVec.z *= z;
		double multByTranspose = recipVec.x + recipVec.y + recipVec.z;
		aniso = exp((2 * M_PI * M_PI) * -(multByTranspose));

		return aniso * me->_occupancy;
	}

	double distSq =	(x * x + y * y + z * z);

	double bf = me->_bFactor;
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
	double scale = 2 * MAX_SCATTERING_DSTAR;

	prepareDistribution(n, scale, this, Absolute::getExpValue);

	return getDistributionCopy();
}

std::vector<BondSample> *Absolute::getManyPositions(BondSampleStyle style)
{
	std::vector<BondSample> *bondSamples = &_bondSamples;
	bondSamples->clear();

	/* B factor isotropic only atm, get mean square displacement in
	 * each dimension. */
	double meanSqDisp = getBFactor() / (8 * M_PI * M_PI);
	meanSqDisp = pow(meanSqDisp, 1./2.);
	double total = 25;

	for (int i = 0; i < total; i++)
	{
		double x = random_norm_dist(0, meanSqDisp);
		double y = random_norm_dist(0, meanSqDisp);
		double z = random_norm_dist(0, meanSqDisp);

		vec3 xyz = make_vec3(x, y, z);
		vec3 full = vec3_add_vec3(xyz, _position);
		BondSample sample;
		sample.basis = make_mat3x3();
		sample.occupancy = 1 / total;
		sample.torsion = 0;
		sample.old_start = make_vec3(0, 0, 0);
		sample.start = full;
		bondSamples->push_back(sample);
	}

	return bondSamples;
}

void Absolute::addToMonomer(MonomerPtr monomer)
{
	makeAtom();
	monomer->addAtom(_atom);
	monomer->getPolymer()->addAtom(_atom);

	Model::addToMonomer(monomer);
}

vec3 Absolute::getStaticPosition()
{
	return _position;
}
