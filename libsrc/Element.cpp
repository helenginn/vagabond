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

std::vector<ElementPtr> Element::elements;

void Element::setupElements()
{
	elements.push_back(ElementPtr(new Element("C", "carbon", 0.70, 6)));
	elements.push_back(ElementPtr(new Element("N", "nitrogen", 0.65, 7)));
	elements.push_back(ElementPtr(new Element("O", "oxygen", 0.60, 8)));
	elements.push_back(ElementPtr(new Element("S", "sulphur", 1.00, 16)));
}

Element::Element(std::string symbol, std::string name, double covRadius,
				 double mass)
{
	_symbol = symbol;
	_name = name;
	_covRadius = covRadius;
	_mass = mass;
	_density = _mass / pow(_covRadius, 3);
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

	double scale = ATOM_SAMPLING;
	double radius = ATOM_MAX_RADIUS; // Angs^-1.

	_shape = FFTPtr(new cFFTW3d());
	_shape->create(2 * radius / scale);
	_shape->setScales(scale);

	double radiusSq = _covRadius * _covRadius;

	for (double x = -radius; x <= radius; x += scale)
	{
		double xfrac = x / (2 * radius);

		for (double y = -radius; y <= radius; y += scale)
		{
			double yfrac = y / (2 * radius);

			for (double z = -radius; z <= radius; z += scale)
			{
				double zfrac = z / (2 * radius);

				double distSq = x * x + y * y + z * z;

				double val = 0;

				if (distSq < radiusSq)
				{
					val = _density / 20 * (1 - distSq);
				}

				_shape->setReal(xfrac, yfrac, zfrac, val);
			}
		}
	}

	_shape->createFFTWplan(1, false);
	_shape->fft(1);
	return _shape;
}