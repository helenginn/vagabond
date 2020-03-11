// Vagabond
// Copyright (C) 2017-2018 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#include "Absolute.h"
#include "shared_ptrs.h"
#include "Atom.h"
#include <math.h>
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
	initialise();
}

void Absolute::initialise()
{
	_position = make_vec3(0, 0, 0);
	_bFactor = 0;
	_element = "";
	_occupancy = 1;
	_hetatm = false;
	_usingTensor = false;
	_tensor = make_mat3x3();
}

mat3x3 Absolute::getRealSpaceTensor()
{
	if (!_usingTensor)
	{
		mat3x3 realSpaceTensor = make_mat3x3();
		double scale = b2var(getMeanSquareDeviation());
		mat3x3_mult_scalar(&realSpaceTensor, scale);
		return realSpaceTensor;
	}
	
	double subtract = 0;
	double mult = 1;
	if (getMolecule())
	{
		subtract = getMolecule()->getAbsoluteBFacSubt();
		mult = getMolecule()->getAbsoluteBFacMult();
	}

	subtract = b2var(subtract);
	mat3x3 copy = _realSpaceTensor;
	mat3x3_mult_scalar(&copy, mult);
	copy.vals[0] -= subtract;
	copy.vals[4] -= subtract;
	copy.vals[8] -= subtract;

	return copy;
}

Absolute::Absolute(vec3 pos, double bFac, 
                   std::string element, double occValue)
{
	initialise();
	_position = pos;
	_bFactor = bFac;
	_element = element;
	trim(_element);
	_occupancy = occValue;
	_hetatm = false;
	_usingTensor = false;
	_tensor = make_mat3x3();
}

std::string Absolute::shortDesc()
{
	std::ostringstream ss;
	ss << "Abs_" + getAtom()->shortDesc();
	return ss.str();
}

AtomPtr Absolute::makeAtom()
{
	AtomPtr myAtom = AtomPtr(new Atom());
	myAtom->setModel(shared_from_this());
	myAtom->setInitialPosition(_position);
	myAtom->setInitialBFactor(_bFactor);
	myAtom->setPDBPosition(_position);
	myAtom->setHetatm(_hetatm);
	myAtom->setAtomNum(_atomNum);
	myAtom->setAlternativeConformer(_conformer);
	ElementPtr element = Element::getElement(_element);
	
	myAtom->setElement(element);
	myAtom->setAtomName(_atomName);
	myAtom->findAtomType(_resName);
	myAtom->setOriginalOccupancy(_occupancy);

	_atom = myAtom;
	return myAtom;
}

void Absolute::addToMolecule(MoleculePtr molecule)
{
	AtomPtr atom = makeAtom();
	molecule->addAtom(atom);
	setMolecule(molecule);
}

double Absolute::getExpValue(void *object, double x, double y, double z)
{
	Absolute *me = static_cast<Absolute *>(object);
	double aniso = 0;
	double mult = 1;
	MoleculePtr molecule = me->getMolecule();
	double subtract = 0;

	if (me->hasMolecule())
	{
		mult = molecule->getAbsoluteBFacMult();
		subtract = molecule->getAbsoluteBFacSubt();
	}

	double sub = b2var(subtract);

	if (me->_usingTensor)
	{
		mat3x3 scaledTensor = me->getRealSpaceTensor();
		vec3 recipVec = make_vec3(x, y, z);
		mat3x3_mult_vec(scaledTensor, &recipVec);

		recipVec.x *= x;
		recipVec.y *= y;
		recipVec.z *= z;
		double multByTranspose = recipVec.x + recipVec.y + recipVec.z;

		// it is 2 * M_PI * M_PI, not 8.
		aniso = exp((2 * M_PI * M_PI) * -(multByTranspose));

		return aniso;
	}

	double distSq = (x * x + y * y + z * z);

	double bf = me->getMeanSquareDeviation();

	if (me->hasMolecule())
	{
		double subtract = molecule->getAbsoluteBFacSubt();
	}

	double exponent = (-0.25) * bf * distSq;
	double value = exp(exponent);

	return value;
}

void Absolute::resetSamples()
{
	_bondSamples.clear();
	_recalcFinal = true;
}

void Absolute::addToMonomer(MonomerPtr monomer)
{
	AtomPtr newAtom = makeAtom();

	monomer->addAtom(newAtom);
	setMolecule(monomer->getPolymer());

	Model::addToMonomer(monomer);
}

double Absolute::getMeanSquareDeviation()
{
	double b = _bFactor;

	if (hasMolecule())
	{
		double subtract = getMolecule()->getAbsoluteBFacSubt();
		double mult = getMolecule()->getAbsoluteBFacMult();
		
		b = mult * (b - subtract);
	}

	return b;
}

void Absolute::setTensor(mat3x3 tensor)
{
	_tensor = tensor;
	_usingTensor = true;

	_realSpaceTensor = _tensor;

	Anisotropicator tropicator;
	tropicator.setTensor(_realSpaceTensor);
	vec3 longestAxis = tropicator.longestAxis();

	getAtom()->setEllipsoidLongestAxis(longestAxis);
	getAtom()->setTensor(_tensor);
}

void Absolute::addProperties()
{
	addDoubleProperty("occupancy", &_occupancy);
	addStringProperty("chain", &_chainID);
	addStringProperty("element", &_element);
	addStringProperty("res_name", &_resName);
	addStringProperty("atom_name", &_atomName);
	addStringProperty("conformer", &_conformer);
	addIntProperty("res_num", &_resNum);
	addIntProperty("atom_num", &_atomNum);
	addBoolProperty("hetatm", &_hetatm);
	addBoolProperty("using_tensor", &_usingTensor);
	addVec3Property("position", &_position);
	addMat3x3Property("tensor", &_tensor);
	addDoubleProperty("bfactor", &_bFactor);

	addReference("atom", _atom.lock());

	Model::addProperties();
}

void Absolute::postParseTidy()
{
	_realSpaceTensor = _tensor;
}

void Absolute::linkReference(BaseParserPtr object, std::string category)
{
	if (category == "atom")
	{
		AtomPtr atom = ToAtomPtr(object);
		_atom = atom;
	}
}

vec3 Absolute::getRandomPosition()
{
	/** Assuming isotropic */
	vec3 randvec = make_vec3(0, 0, 0);	
	double stdev = sqrt(b2var(_bFactor));

	for (int i = 0; i < 3; i++)
	{
		double deviation = random_norm_dist(0, stdev);
		*(&randvec.x + i) = deviation;
	}

	vec3 absPos = getAbsolutePosition();
	vec3 total = vec3_add_vec3(absPos, randvec);

	return total;
}

