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

Atom::Atom()
{

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

void Atom::addToMap(FFTPtr fft, mat3x3 unit_cell)
{
	FFTPtr atomDist = _element->getDistribution();
	FFTPtr modelDist = getBlur();
	FFTPtr modified = std::make_shared<cFFTW3d>(*modelDist);

	cFFTW3d::multiply(modified, atomDist);
	modified->fft(1);

	double xPos = getPosition().x;
	double yPos = getPosition().y;
	double zPos = getPosition().z;

	vec3 pos = make_vec3(xPos, yPos, zPos);
	mat3x3_mult_vec(unit_cell, &pos);

	cFFTW3d::add(fft, modified, 2, pos.x, pos.y, pos.z, false, MaskProtein);
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