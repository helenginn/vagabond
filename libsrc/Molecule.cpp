//
//  Molecule.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "shared_ptrs.h"
#include "Molecule.h"
#include "Atom.h"
#include "Element.h"
#include "Bond.h"
#include <float.h>
#include <iostream>
#include "fftw3d.h"
#include "mat3x3.h"

void Molecule::addModel(ModelPtr aModel)
{
	models.push_back(aModel);
}

void Molecule::addAtom(AtomPtr atom)
{
	atoms.push_back(atom);
}

void Molecule::addToMap(FFTPtr fft, mat3x3 _real2frac)
{
	for (int i = 0; i < atomCount(); i++)
	{
		atoms[i]->addToMap(fft, _real2frac);
	}
}

void Molecule::summary()
{
	std::cout << "| I am chain " << getChainID() << std::endl;
	std::cout << "| Atoms: " << atomCount() << std::endl;
}

void Molecule::refine(CrystalPtr target, RefinementType rType)
{

}

void Molecule::makePDB(std::string filename)
{

}

void Molecule::tiedUpScattering(double *tied, double *all)
{
	double total = 0;
	double some = 0;
	double someCount = 0;
	double someTotal = 0;

	for (int i = 0; i < atomCount(); i++)
	{
		if (atoms[i]->getModel()->getClassName() == "Bond")
		{
			some += atoms[i]->getElement()->electronCount();
			someCount++;
		}

		total += atoms[i]->getElement()->electronCount();
		someTotal++;
	}

	double scatterFrac = sqrt(some / total);

	*tied += some;
	*all += total;
}

void Molecule::resetInitialPositions()
{
	for (int i = 0; i < atomCount(); i++)
	{
		atom(i)->setInitialPosition(atom(i)->getPosition());

		ModelPtr model = atom(i)->getModel();

		if (model->isBond())
		{
			ToBondPtr(model)->resetBondDirection();
		}
	}
}

