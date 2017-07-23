//
//  Element.cpp
//  vagabond
//
//  Created by Helen Ginn on 16/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Element.h"
#include <math.h>
#include "Shouter.h"
#include "fftw3d.h"
#include <stdlib.h>
#include <iostream>

std::vector<ElementPtr> Element::elements;

void Element::setupElements()
{
	elements.push_back(ElementPtr(new Element("H", "hydrogen")));
	elements.push_back(ElementPtr(new Element("C", "carbon")));
	elements.push_back(ElementPtr(new Element("N", "nitrogen")));
	elements.push_back(ElementPtr(new Element("O", "oxygen")));
	elements.push_back(ElementPtr(new Element("S", "sulphur")));
}

Element::Element(std::string symbol, std::string name)
{
	_symbol = symbol;
	_name = name;

	if (_symbol == "H")
	{
		memcpy(_scattering, ScatterFactors::hScatter, ScatterFactors::numScatter * sizeof(float));
	}
	else if (_symbol == "C")
	{
		memcpy(_scattering, ScatterFactors::cScatter, ScatterFactors::numScatter * sizeof(float));
	}
	else if (_symbol == "N")
	{
		memcpy(_scattering, ScatterFactors::nScatter, ScatterFactors::numScatter * sizeof(float));
	}
	else if (symbol == "O")
	{
		memcpy(_scattering, ScatterFactors::oScatter, ScatterFactors::numScatter * sizeof(float));
	}
	else if (_symbol == "S")
	{
		memcpy(_scattering, ScatterFactors::sScatter, ScatterFactors::numScatter * sizeof(float));
	}
	else
	{
		memset(_scattering, 0, ScatterFactors::numScatter * sizeof(float));
	}
}

ElementPtr Element::getElement(std::string symbol)
{
	if (!elements.size())
	{
		setupElements();
	}

	for (int i = 0; i < elements.size(); i++)
	{
		if (elements[i]->getSymbol() == symbol)
		{
			return elements[i];
		}
	}

	shout_at_user("Missing element of symbol " + symbol);

	return ElementPtr();
}


/*  Problems: value needs to be proportional to the density of the atom
 *  Need to do a better calculation of actual voxels.
 */
FFTPtr Element::getDistribution()
{
	if (_shape)
	{
		return _shape;
	}

	_shape = FFTPtr(new cFFTW3d());
	double n = ATOM_SAMPLING_COUNT / PROTEIN_SAMPLING;
	double scale = MAX_SCATTERING_DSTAR / n;

	_shape->create(n);
	_shape->setScales(scale * PROTEIN_SAMPLING);
//	double sampling = 1 / (scale * n);

	if (_scattering[0] <= 0)
	{
		shout_at_user("Helen needs to add your element\n" \
					  + _symbol + " to the list of atoms.\n" \
					  "Sorry!\n");
	}

	int totalScatterPoints = ScatterFactors::numScatter;
	
	for (double x = -0.5; x <= 0.5; x += 1 / n)
	{
		for (double y = -0.5; y <= 0.5; y += 1 / n)
		{
			for (double z = -0.5; z <= 0.5; z += 1 / n)
			{
				double mod = MAX_SCATTERING_DSTAR * 2;
				double xAng = x * mod;
				double yAng = y * mod;
				double zAng = z * mod;

				double distSq =	(xAng * xAng + yAng * yAng + zAng * zAng);
				double dist = sqrt(distSq);

				double val = 0;

				for (int i = 0; i < totalScatterPoints - 1; i++)
				{
					if ((dist >= ScatterFactors::dScatter[i] &&
						dist < ScatterFactors::dScatter[i + 1]))
					{
						double interpolateToNext = (ScatterFactors::dScatter[i] - dist);
						interpolateToNext /= fabs(ScatterFactors::dScatter[i + 1] - ScatterFactors::dScatter[i]);
						val = (_scattering[i] + interpolateToNext * (_scattering[i] - _scattering[i + 1]));
						break;
					}
				}

				_shape->setReal(x, y, z, val);
			}
		}
	}

	_shape->createFFTWplan(1, false);

	return _shape;
}