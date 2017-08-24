//
//  Atom.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Atom.h"
#include "fftw3d.h"
#include "mat3x3.h"
#include <math.h>
#include <stdlib.h>
#include "Model.h"
#include <iostream>
#include "Element.h"
#include "Monomer.h"
#include "Polymer.h"
#include <iomanip>
#include "FileReader.h"
#include <sstream>

Atom::Atom()
{
	_initialPosition = make_vec3(0, 0, 0);
	_initialB = 0;
	_geomType = AtomUnassigned;
}

Atom::Atom(Atom &other)
{
	_initialPosition = other._initialPosition;
	_initialB = other._initialB;
	_geomType = other._geomType;
	_element = other._element;
	_atomName = other._atomName;
	_model = other._model;
	_monomer = other._monomer;
}

void Atom::inheritParents()
{
	getMonomer()->addAtom(shared_from_this());
	getMonomer()->getPolymer()->addAtom(shared_from_this());
}

void Atom::setModel(ModelPtr model)
{
	_model = model;
}

FFTPtr Atom::getBlur()
{
	FFTPtr modelDist = _model->getDistribution();

	if (_distModelOnly)
	{
		modelDist = _distModelOnly->getDistribution();
	}

	return modelDist;
}

double Atom::scoreWithMap(FFTPtr fft, mat3x3 unit_cell,
						  std::vector<double> *xs, std::vector<double> *ys)
{
	FFTPtr atomDist = _element->getDistribution();
	FFTPtr modelDist = getBlur();
	FFT::multiply(modelDist, atomDist);
	modelDist->fft(1);
	modelDist->invertScale();

	vec3 pos = getPosition();
	mat3x3_mult_vec(unit_cell, &pos);

	double score = FFT::score(fft, modelDist, pos, xs, ys);

	return score;
}

void Atom::addToMap(FFTPtr fft, mat3x3 unit_cell, vec3 offset)
{
	FFTPtr atomDist = _element->getDistribution();
	FFTPtr modified = getBlur();

	FFT::multiply(modified, atomDist);
	modified->fft(1);
//	modified->printSlice();
	modified->invertScale();

	vec3 pos = _model->getAbsolutePosition();
	pos = vec3_subtract_vec3(pos, offset);
	mat3x3_mult_vec(unit_cell, &pos);

	FFT::add(fft, modified, pos);
}

vec3 Atom::getPosition()
{
	return _model->getStaticPosition();
}

bool Atom::isBackbone()
{
	if (_atomName == "C") return true;
	if (_atomName == "H") return true;
	if (_atomName == "N") return true;
	if (_atomName == "O") return true;
	if (_atomName == "HA3") return true;

	return false;
}

bool Atom::isBackboneAndSidechain()
{
	if (_atomName == "CA") return true;
	if (_atomName == "HA") return true;
	if (_atomName == "HA2") return true;

	return false;
}

/* Convert to lookup table */
void Atom::findAtomType(std::string resName)
{
	_geomType = GeomTable::getGeomTable().getType(resName, _atomName);
}

std::string Atom::pdbLineBeginning(int i)
{
	std::string residueName = getMonomer()->getIdentifier();
	int resNum = getMonomer()->getResidueNum();
	to_upper(residueName);
	std::ostringstream line;

	char conformer[] = "A";
	conformer[0] += i;

	line << "ATOM  ";
	line << std::setfill(' ') << std::setw(5) << std::fixed << _atomNum;
	line << std::setfill(' ') << std::setw(4) << _atomName;
	line << " " << conformer;
	line << std::setw(3) << residueName;
	line << " A";
	line << std::setfill(' ') << std::setw(4) << resNum;
	line << "    ";

	return line.str();
}

void Atom::setKeepModel()
{
	_distModelOnly = _model;
}