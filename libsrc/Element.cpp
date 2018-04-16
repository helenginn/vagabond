//
//  Element.cpp
//  vagabond
//
//  Created by Helen Ginn on 16/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Element.h"
#include "Atom.h"
#include <math.h>
#include "Shouter.h"
#include "fftw3d.h"
#include <stdlib.h>
#include <iostream>
#include "../libinfo/ScatterFactors.h"
#include "Options.h"

std::vector<ElementPtr> Element::elements;

void Element::setupElements()
{
	elements.push_back(ElementPtr(new Element("H", "hydrogen", 1, ScatterFactors::hScatter)));
	elements.push_back(ElementPtr(new Element("C", "carbon", 6,  ScatterFactors::cScatter)));
	elements.push_back(ElementPtr(new Element("N", "nitrogen", 7,  ScatterFactors::nScatter)));
	elements.push_back(ElementPtr(new Element("O", "oxygen", 8,  ScatterFactors::oScatter)));
	elements.push_back(ElementPtr(new Element("F", "fluorine", 9,  ScatterFactors::fScatter)));
	elements.push_back(ElementPtr(new Element("NA", "sodium", 11,  ScatterFactors::naScatter)));
	elements.push_back(ElementPtr(new Element("MG", "magnesium", 12,  ScatterFactors::mgScatter)));
	elements.push_back(ElementPtr(new Element("P", "phosphorus", 15,  ScatterFactors::pScatter)));
	elements.push_back(ElementPtr(new Element("S", "sulphur", 16,  ScatterFactors::sScatter)));
	elements.push_back(ElementPtr(new Element("CL", "chlorine", 17,  ScatterFactors::clScatter)));
	elements.push_back(ElementPtr(new Element("CA", "calcium", 20,  ScatterFactors::caScatter)));
	elements.push_back(ElementPtr(new Element("MN", "manganese", 25,  ScatterFactors::mnScatter)));
	elements.push_back(ElementPtr(new Element("FE", "iron", 26,  ScatterFactors::feScatter)));
	elements.push_back(ElementPtr(new Element("ZN", "zinc", 30,  ScatterFactors::znScatter)));
}

Element::Element(std::string symbol, std::string name, double electrons, const float *scatter)
{
	_symbol = symbol;
	_name = name;
	_electrons = electrons;
	memcpy(_scattering, scatter, ScatterFactors::numScatter * sizeof(float));

	if (_symbol == "S")
	{
		for (int i = 0; i < ScatterFactors::numScatter; i++)
		{
			_scattering[i] += 0.2;
		}
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

	shout_at_user("Missing element of symbol \"" + symbol + "\"");

	return ElementPtr();
}

double Element::getSolventMaskValue(void *object, double x, double y, double z)
{
	Element *me = static_cast<Element *>(object);
	
	if (me->electronCount() <= 1)
	{
		return 0;
	}

	double distSq = (x * x + y * y + z * z);
	double dist = sqrt(distSq);
	
	/** Point at which solvent should be completely disallowed */
	double min = 1.8;
	/** Point at which solvent is completely free to exist */
	double max = 2.5;
	
	if (dist <= min)
	{
		return 1;
	}
	else if (dist > max)
	{
		return 0;	
	}
	else
	{
		double extra = dist - min;
		extra /= (max - min);
		
		return 1 - extra;
	}
}

double Element::getVoxelValue(void *object, double x, double y, double z)
{
	Element *me = static_cast<Element *>(object);
	int totalScatterPoints = ScatterFactors::numScatter;

	double distSq = (x * x + y * y + z * z);
	double dist = sqrt(distSq);
	dist *= Options::getProteinSampling();

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

FFTPtr Element::getMask()
{
	get_voxel_value *func = getSolventMaskValue;	
	double n = 12;
	double scale = Options::getProteinSampling();

	FFTPtr fft = FFTPtr(new FFT());
	fft->create(n);
	fft->setScales(scale);

	for (int x = -fft->nx / 2; x <= fft->nx / 2; x++)
	{
		for (int y = -fft->ny / 2; y <= fft->ny / 2; y++)
		{
			for (int z = -fft->nz / 2; z <= fft->nz / 2; z++)
			{
				double xAng = x * scale;
				double yAng = y * scale;
				double zAng = z * scale;

				double val = getSolventMaskValue(this, xAng, yAng, zAng);
				long element = fft->element(x, y, z);
				fft->setElement(element, val, 0);
			}
		}
	}

	return fft;
}

FFTPtr Element::getDistribution(bool, int new_n)
{
	double n = ATOM_SAMPLING_COUNT;
	
	if (new_n > 0)
	{
		n = new_n;	
	}
	
	double maxDStar = Options::getRuntimeOptions()->getActiveCrystalDStar();
	double scale = 2.0 * maxDStar;
	
	get_voxel_value *func = getVoxelValue;

	prepareDistribution(n, scale, this, func);
	return getDistributionCopy();
}


std::vector<ElementPtr> Element::elementList(std::vector<AtomPtr> atoms)
{
	std::vector<ElementPtr> elements;
	
	for (size_t i = 0; i < atoms.size(); i++)
	{
		ElementPtr element = atoms[i]->getElement();
		std::vector<ElementPtr>::iterator it;
		it = std::find(elements.begin(), elements.end(), element);
		
		if (it == elements.end())
		{
			elements.push_back(element);
		}
	}
	
	return elements;
}

