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

void Atom::addToMap(FFTPtr fft, mat3x3 unit_cell)
{
	FFTPtr atomDist = _element->getDistribution();
	FFTPtr modified = getBlur();

	FFT::multiply(modified, atomDist);
	modified->fft(1);
	modified->invertScale();

	vec3 pos = getPosition();
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
	if (_atomName == "H") return true;

	return false;
}

bool Atom::isBackboneAndSidechain()
{
	if (_atomName == "CA") return true;

	return false;
}

/* Convert to lookup table */
void Atom::findAtomType(std::string resName)
{
	if (_atomName == "CA" && resName != "gly")
	{
		_geomType = AtomCH1E;
	}
	else if (_atomName == "CB" && resName == "thr")
	{
		_geomType = AtomCH1E;
	}
	else if (_atomName == "CG2" && resName == "thr")
	{
		_geomType = AtomCH3E;
	}
	else if (_atomName == "OG1" && resName == "thr")
	{
		_geomType = AtomOH1;
	}
	else if (_element->getSymbol() == "C" && resName == "lys")
	{
		_geomType = AtomCH2E;
	}
	else if (_element->getSymbol() == "N" && resName == "lys")
	{
		_geomType = AtomNH3;
	}
}