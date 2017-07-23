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

void Atom::addConnection(ModelPtr model)
{
	connections.push_back(model);
}


FFTPtr Atom::addToMap(FFTPtr fft, FFTPtr reuseModelDist, mat3x3 unit_cell)
{
	FFTPtr atomDist = _element->getDistribution();

	connections[0]->getDistribution(&reuseModelDist);
	cFFTW3d::multiply(reuseModelDist, atomDist);
	reuseModelDist->fft(1);

	double xPos = getPosition().x;
	double yPos = getPosition().y;
	double zPos = getPosition().z;

	vec3 pos = make_vec3(xPos, yPos, zPos);
	mat3x3_mult_vec(unit_cell, &pos);

	cFFTW3d::add(fft, reuseModelDist, 2, pos.x, pos.y, pos.z, false, MaskProtein);

	return reuseModelDist;
}