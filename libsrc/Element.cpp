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
#include "ScatterFactors.h"

std::vector<ElementPtr> Element::elements;

void Element::setupElements()
{
	elements.push_back(ElementPtr(new Element("H", "hydrogen", 1)));
	elements.push_back(ElementPtr(new Element("C", "carbon", 6)));
	elements.push_back(ElementPtr(new Element("N", "nitrogen", 7)));
	elements.push_back(ElementPtr(new Element("O", "oxygen", 8)));
	elements.push_back(ElementPtr(new Element("S", "sulphur", 16)));
	elements.push_back(ElementPtr(new Element("CL", "chlorine", 17)));
}

Element::Element(std::string symbol, std::string name, double electrons)
{
	_symbol = symbol;
	_name = name;
	_electrons = electrons;

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
	else if (_symbol == "CL")
	{
		memcpy(_scattering, ScatterFactors::clScatter, ScatterFactors::numScatter * sizeof(float));
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

double Element::getVoxelValue(void *object, double x, double y, double z)
{
	Element *me = static_cast<Element *>(object);
	int totalScatterPoints = ScatterFactors::numScatter;

	double distSq =	(x * x + y * y + z * z);
	double dist = sqrt(distSq);
	dist *= PROTEIN_SAMPLING;

	double val = 0;

	for (int i = 0; i < totalScatterPoints - 1; i++)
	{
		if ((dist >= ScatterFactors::dScatter[i] &&
			 dist < ScatterFactors::dScatter[i + 1]))
		{
			double interpolateToNext = (ScatterFactors::dScatter[i] - dist);
			interpolateToNext /= fabs(ScatterFactors::dScatter[i + 1] - ScatterFactors::dScatter[i]);
			val = (me->_scattering[i] + interpolateToNext * (me->_scattering[i] - me->_scattering[i + 1]));
			break;
		}
	}

	return val;
}

FFTPtr Element::getDistribution()
{
	double n = ATOM_SAMPLING_COUNT;
	double scale = 2 * MAX_SCATTERING_DSTAR;

	return prepareDistribution(n, scale, this, getVoxelValue);
}