//
//  Absolute.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Absolute.h"
#include "shared_ptrs.h"
#include "Molecule.h"
#include "Atom.h"
#include <math.h>
#include "fftw3d.h"
#include "Element.h"
#include <iostream>

Absolute::Absolute(vec3 pos, double bFac, std::string element, double occValue)
{
	position = pos;
	bFactor = bFac;
	_element = element;
	_occupancy = occValue;
}

void Absolute::addToMolecule(MoleculePtr molecule)
{
	Model::addToMolecule(molecule);

	AtomPtr myAtom = AtomPtr(new Atom());
	myAtom->addConnection(shared_from_this());
	ElementPtr element = Element::getElement(_element);
	myAtom->setElement(element);

	myAtom->setPosition(position);
	_atom = myAtom;

	molecule->addAtom(myAtom);
}

/*  Absolute distribution only needs to be the blurring due to the atomic
 *  B factor. The position should be provided by a different function.
 */
FFTPtr Absolute::getDistribution()
{
	double scale = ATOM_SAMPLING;
	double radius = ATOM_MAX_RADIUS; // Angs^-1.

	FFTPtr fft = FFTPtr(new cFFTW3d());
	int n = 2 * radius / scale;
	fft->create(n);
	fft->setScales(scale);
	double sampling = 1 / (ATOM_SAMPLING * n);
	double switch_sampling = sampling / scale;

	for (double x = -radius; x <= radius; x += scale)
	{
		double xfrac = x / (2 * radius);

		for (double y = -radius; y <= radius; y += scale)
		{
			double yfrac = y / (2 * radius);

			for (double z = -radius; z <= radius; z += scale)
			{
				double zfrac = z / (2 * radius);

				double xAng = x + scale / 2;
				double yAng = y + scale / 2;
				double zAng = z + scale / 2;

				double distSq = switch_sampling * switch_sampling *
			//	(xAng * xAng + yAng * yAng + zAng * zAng);
				x * x + y * y + z * z;

				double value = exp((-0.25) * bFactor * distSq);
				value *= _occupancy;

				fft->setReal(xfrac, yfrac, zfrac, value);
			}
		}
	}

	fft->createFFTWplan(1, false);
	
	return fft;
}