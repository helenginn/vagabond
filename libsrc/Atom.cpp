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

double Atom::scoreWithMap(FFTPtr fft, mat3x3 unit_cell)
{
	FFTPtr atomDist = _element->getDistribution();
	FFTPtr modelDist = getBlur();
	FFT::multiply(modelDist, atomDist);
	modelDist->fft(1);

	vec3 pos = getPosition();
	mat3x3_mult_vec(unit_cell, &pos);

	double score = FFT::score(fft, modelDist, pos);

	return score;
}

void Atom::addToMap(FFTPtr fft, mat3x3 unit_cell)
{
	FFTPtr atomDist = _element->getDistribution();
	FFTPtr modelDist = getBlur();
	FFTPtr modified = std::make_shared<FFT>(*modelDist);

	FFT::multiply(modified, atomDist);
	modified->fft(1);

	vec3 pos = getPosition();
	mat3x3_mult_vec(unit_cell, &pos);

	FFT::add(fft, modified, 2, pos.x, pos.y, pos.z, false, MaskProtein);
}

vec3 Atom::getPosition()
{
	return _model->getPosition();
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