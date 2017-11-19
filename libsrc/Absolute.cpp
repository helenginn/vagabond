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
#include "Crystal.h"
#include "Anisotropicator.h"
#include <iomanip>
#include <iostream>

Absolute::Absolute()
{
	_position = make_vec3(0, 0, 0);
	_bFactor = 0;
	_element = "";
	_occupancy = 1;
	_hetatm = false;
	_usingTensor = false;
	_tensor = make_mat3x3();
}

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
	myAtom->setPDBPosition(_position);
	myAtom->setAtomNum(_atomNum);
	myAtom->setAlternativeConformer(_conformer);
	ElementPtr element = Element::getElement(_element);
	myAtom->setElement(element);
	myAtom->setAtomName(_atomName);
	myAtom->findAtomType(_resName);
	myAtom->setOriginalOccupancy(_occupancy);

	_atom = myAtom;
}

void Absolute::addToMolecule(MoleculePtr molecule)
{
	makeAtom();

	molecule->addAtom(_atom);

	setMolecule(molecule);
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

	if (me->hasMolecule())
	{
		MoleculePtr molecule = me->getMolecule();
		double subtract = molecule->getAbsoluteBFacSubtract();
		double mult = molecule->getAbsoluteBFacMult();
		bf -= subtract;
		bf *= mult;
	}

	double exponent = (-0.25) * bf * distSq;
	double value = exp(exponent);
	value *= me->_occupancy;

	return value;
}

/*  Absolute distribution only needs to be the blurring due to the atomic
 *  B factor. The position should be provided by a different function. */

FFTPtr Absolute::getDistribution(bool quick)
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

	if (style == BondSampleStatic)
	{
		BondSample sample;
		sample.basis = make_mat3x3();
		sample.occupancy = 1;
		sample.torsion = 0;
		sample.old_start = make_vec3(0, 0, 0);
		sample.start = _position;
		bondSamples->push_back(sample);
		return bondSamples;
	}

		/* B factor isotropic only atm, get mean square displacement in
	 * each dimension. */
	double meanSqDisp = getBFactor() / (8 * M_PI * M_PI);
	meanSqDisp = sqrt(meanSqDisp);

	double occTotal = 0;

	int samples = 81;
	int rnd = 1;
	double total = 2;
	double step = (meanSqDisp * 1.5) / total; // cover four st.dev.s

	std::vector<vec3> points;
	double offset = 2. / (double)samples;
	double increment = M_PI * (3.0 - sqrt(5));

	for (int i = 0; i < samples; i++)
	{
		double y = (((double)i * offset) - 1) + (offset / 2);
		double r = sqrt(1 - y * y);
		r *= meanSqDisp;

		double phi = (double)((i + rnd) % samples) * increment;

		double x = cos(phi) * r;
		double z = sin(phi) * r;
		y *= r;

		points.push_back(make_vec3(x, y, z));
	}

	for (int i = 0; i < points.size(); i++)
	{
		vec3 full = vec3_add_vec3(points[i], _position);
		double occ = 1;
		occTotal += occ;

		BondSample sample;
		sample.basis = make_mat3x3();
		sample.occupancy = occ;
		sample.torsion = 0;
		sample.old_start = make_vec3(0, 0, 0);
		sample.start = full;
		bondSamples->push_back(sample);
	}

/*

    for (int i = -total; i <= total; i++)
	{
		double xVal = i * step;

		for (int j = -total; j <= total; j++)
		{
			double yVal = j * step;

			for (int k = -total; k <= total; k++)
			{
				double zVal = k * step;

				double distance = sqrt(xVal * xVal + yVal * yVal + zVal * zVal);
				double stdev = distance / meanSqDisp;

				if (stdev > 2.00) continue;

				double occ = normal_distribution(xVal, meanSqDisp);
				occ *= normal_distribution(yVal, meanSqDisp);
				occ *= normal_distribution(zVal, meanSqDisp);
				occ = 1;
				occTotal += occ;
				vec3 xyz = make_vec3(xVal, yVal, zVal);
				vec3 full = vec3_add_vec3(xyz, _position);

				BondSample sample;
				sample.basis = make_mat3x3();
				sample.occupancy = occ;
				sample.torsion = 0;
				sample.old_start = make_vec3(0, 0, 0);
				sample.start = full;
				bondSamples->push_back(sample);
			}
		}
	}
*/

	for (int i = 0; i < bondSamples->size(); i++)
	{
		(*bondSamples)[i].occupancy /= occTotal;
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

double Absolute::getMeanSquareDeviation()
{
	return _bFactor;
}

vec3 Absolute::getStaticPosition()
{
	return _position;
}

void Absolute::setTensor(mat3x3 tensor, CrystalPtr crystal)
{
	_tensor = tensor;
	_usingTensor = true;

	mat3x3 toNormal = crystal->getHKL2Real();
	_realSpaceTensor = mat3x3_mult_mat3x3(toNormal, _tensor);

	Anisotropicator tropicator;
	tropicator.setTensor(_realSpaceTensor);
	vec3 longestAxis = tropicator.longestAxis();

	getAtom()->setEllipsoidLongestAxis(longestAxis);
}
