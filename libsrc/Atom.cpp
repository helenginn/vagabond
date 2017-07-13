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
	double sampling = fft->getSampling();
	
	double xPos = getPosition().x * sampling;
	double yPos = getPosition().y * sampling;
	double zPos = getPosition().z * sampling;

	std::cout << "pos = " << xPos << ", " << yPos << ", " << zPos << std::endl;

	for (int k = -1; k < 2; k++)
	{
		for (int j = -1; j < 2; j++)
		{
			for (int i = -1; i < 2; i++)
			{
				vec3 pos = make_vec3(xPos + (double)i,
									 yPos + (double)j,
									 zPos + (double)k);
				mat3x3_mult_vec(unit_cell, &pos);

				std::cout << "pos after = " << pos.x << ", " << pos.y << ", " << pos.z << std::endl;

				double val = 4 - abs(i) - abs(j) - abs(k);
				fft->setReal(pos.x, pos.y, pos.z, val);
			}
		}
	}
}