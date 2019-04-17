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
	elements.push_back(ElementPtr(new Element("K", "potassium", 19,  ScatterFactors::kScatter)));
	elements.push_back(ElementPtr(new Element("CA", "calcium", 20,  ScatterFactors::caScatter)));
	elements.push_back(ElementPtr(new Element("MN", "manganese", 25,  ScatterFactors::mnScatter)));
	elements.push_back(ElementPtr(new Element("FE", "iron", 26,  ScatterFactors::feScatter)));
	elements.push_back(ElementPtr(new Element("NI", "nickel", 28,  ScatterFactors::niScatter)));
	elements.push_back(ElementPtr(new Element("CU", "copper", 29,  ScatterFactors::cuScatter)));
	elements.push_back(ElementPtr(new Element("ZN", "zinc", 30,  ScatterFactors::znScatter)));
	elements.push_back(ElementPtr(new Element("SE", "selenium", 34,  ScatterFactors::seScatter)));
	elements.push_back(ElementPtr(new Element("BR", "bromine", 35,  ScatterFactors::brScatter)));
}

Element::Element(std::string symbol, std::string name, double electrons, const float *scatter)
{
	_symbol = symbol;
	_name = name;
	_electrons = electrons;
	memcpy(_scattering, scatter, ScatterFactors::numScatter * sizeof(float));
}

ElementPtr Element::getElement(std::string symbol)
{
	if (!elements.size())
	{
		setupElements();
	}

	for (size_t i = 0; i < elements.size(); i++)
	{
		if (elements[i]->getSymbol() == symbol)
		{
			return elements[i];
		}
	}

	shout_at_user("Missing element of symbol \"" + symbol + "\"");

	return ElementPtr();
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

void Element::populateFFT(mat3x3 basis, FFTPtr fft)
{
	for (int z = 0; z < fft->nz; z++)
	{
		for (int y = 0; y < fft->ny; y++)
		{
			for (int x = 0; x < fft->nx; x++)
			{
				double fx = (double) x / (double)fft->nx;
				double fy = (double) y / (double)fft->ny;
				double fz = (double) z / (double)fft->nz;

				vec3 fr = make_vec3(fx, fy, fz);
				vec3 real = fr;
				mat3x3_mult_vec(basis, &real);
				double val = getVoxelValue(this, real.x, real.y, real.z);
				int ele = fft->elementFromFrac(fr.x, fr.y, fr.z);
				fft->data[ele][0] = val;
			}
		}
	}
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

