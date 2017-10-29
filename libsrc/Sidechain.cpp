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

