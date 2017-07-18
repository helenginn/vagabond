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
	elements.push_back(ElementPtr(new Element("C", "carbon")));
	elements.push_back(ElementPtr(new Element("N", "nitrogen")));
	elements.push_back(ElementPtr(new Element("O", "oxygen")));
	elements.push_back(ElementPtr(new Element("S", "sulphur")));
}

Element::Element(std::string symbol, std::string name)
{
	_symbol = symbol;
	_name = name;

	if (_symbol == "C")
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

	double scale = ATOM_SAMPLING;
	double radius = ATOM_MAX_RADIUS; // Angs^-1.

	_shape = FFTPtr(new cFFTW3d());
	_shape->create(2 * radius / scale);
	_shape->setScales(scale);

	if (_scattering[0] <= 0)
	{
		shout_at_user("Helen needs to add your element\n" \
					  + _symbol + " to the list of atoms.\n" \
					  "Sorry!\n");
	}

	int totalScatterPoints = ScatterFactors::numScatter;

	for (double x = -radius; x <= radius; x += scale)
	{
		double xfrac = x / (2 * radius);

		for (double y = -radius; y <= radius; y += scale)
		{
			double yfrac = y / (2 * radius);

			for (double z = -radius; z <= radius; z += scale)
			{
				double zfrac = z / (2 * radius);

				double xp = x + scale / 2;
				double yp = y + scale / 2;
				double zp = z + scale / 2;

			//	double distSq = xp * xp + yp * yp + zp * zp;
				double distSq = x * x + y * y + z * z;
				double dist = sqrt(distSq);

				double val = 0;

				if (distSq <= 0.005)
				{
					val = _scattering[0];
				}

				for (int i = 0; i < totalScatterPoints - 1; i++)
				{
					if ((dist >= ScatterFactors::dScatter[i] &&
						dist < ScatterFactors::dScatter[i + 1]))
					{
						val = (_scattering[i]);
						break;
					}
				}

				long int element = _shape->elementFromFrac(xfrac, yfrac, zfrac);
				_shape->setElement(element, val, 0);
			}
		}
	}

	_shape->createFFTWplan(1);
	return _shape;
}