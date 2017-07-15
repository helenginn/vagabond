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
#include <iostream>

void Atom::addConnection(ModelPtr model)
{
	connections.push_back(model);
}

void Atom::addToMap(FFTPtr fft, mat3x3 unit_cell)
{
	double xPos = getPosition().x;
	double yPos = getPosition().y;
	double zPos = getPosition().z;

	for (int k = -1; k < 2; k++)
	{
		for (int j = -1; j < 2; j++)
		{
			for (int i = -1; i < 2; i++)
			{
				vec3 pos = make_vec3(xPos + 0.2 * double(i),
									 yPos + 0.2 * double(j),
									 zPos + 0.2 * double(k));
				mat3x3_mult_vec(unit_cell, &pos);

				double val = 4 - abs(i) - abs(j) - abs(k);
				fft->setReal(pos.x, pos.y, pos.z, val);
			}
		}
	}
}