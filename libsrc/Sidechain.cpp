//
//  Sidechain.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Sidechain.h"
#include "Sampler.h"
#include "Bond.h"
#include "Atom.h"
#include <iostream>
#include "Monomer.h"
#include "FileReader.h"
#include "Absolute.h"
#include "Backbone.h"

bool Sidechain::shouldRefineMagicAxis(BondPtr bond)
{
	return (bond->getMinor()->getAtomName() == "CB");
}

void Sidechain::fixBackboneTorsions(AtomPtr betaTorsion)
{
	AtomPtr atom = findAtom("CB");

	if (!atom)
	{
		return;
	}

	ModelPtr model = atom->getModel();

	if (model->isBond())
	{
		ToBondPtr(model)->setTorsionAtoms(betaTorsion);
	}
}

void Sidechain::setInitialDampening()
{
	for (int i = 0; i < atomCount(); i++)
	{
		if (atom(i)->isBackbone())
		{
			continue;
		}

		BondPtr bond = BondPtr();

		if (atom(i)->getModel()->isBond())
		{
			bond = ToBondPtr(atom(i)->getModel());
		}
		else
		{
			continue;
		}

		Bond::setDampening(&*bond, -0.0);

		double kick = 0.6;
		std::string id = getMonomer()->getIdentifier();

		if (id == "tyr" || id == "phe" || id == "trp")
		{
			kick = 0.3;
		}

		if (bond->isRefinable() && atom(i)->getAtomName() == "CB")
		{
			Bond::setTorsionBlur(&*bond, kick);
		}
	}
}

void Sidechain::splitConformers()
{
	int count = conformerCount();

	if (count <= 1) return;

	if (getMonomer()->getBackbone()->findAtoms("N").size() != 1)
	{
		std::cout << "Not splitting whole residue conformer, "
		<< getMonomer()->getResidueNum() << getMonomer()->getIdentifier() << std::endl;
//		return;
	}

	AtomPtr start = findAtom("CB");

	if (!start || !start->getModel()->isBond())
	{
		return;
	}

	BondPtr bond = ToBondPtr(start->getModel());

	if (count > 2) return;

	for (int i = 1; i < count; i++)
	{
		std::string confID = conformer(i);

		bond->splitBond();
	}

	AtomList atoms = findAtoms("CB");

	for (int i = 0; i < atoms.size(); i++)
	{
		AtomPtr atom = atoms[i].lock();
		if (atom->getModel()->isBond())
		{
			BondPtr bond = ToBondPtr(atom->getModel());
			double origOcc = atom->getOriginalOccupancy();
			Bond::setOccupancy(&*bond, origOcc);
		}
	}

	for (int i = 0; i < getMonomer()->atomCount(); i++)
	{
		AtomPtr atom = getMonomer()->atom(i);

		if (atom->getModel()->isAbsolute() && atom->getAlternativeConformer().length())
		{
			atom->setWeighting(0);
		}
	}
}
