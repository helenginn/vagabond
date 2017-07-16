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
	fft->create(2 * radius / scale);
	fft->setScales(scale);

	double B = bFactor * 8 * M_PI;

	for (double x = -radius; x <= radius; x += scale)
	{
		double xfrac = x / (2 * radius);

		for (double y = -radius; y <= radius; y += scale)
		{
			double yfrac = y / (2 * radius);

			for (double z = -radius; z <= radius; z += scale)
			{
				double zfrac = z / (2 * radius);
				double xAng = x * scale;
				double yAng = y * scale;
				double zAng = z * scale;

				double distSq = xAng * xAng + yAng * yAng + zAng * zAng;

				double value = _occupancy * exp((-0.25) * B * distSq);

				fft->setReal(xfrac, yfrac, zfrac, value);
			}
		}
	}

	fft->createFFTWplan(1, false);

	return fft;
}