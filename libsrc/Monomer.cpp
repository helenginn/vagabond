//
//  Monomer.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Monomer.h"
#include "Backbone.h"
#include "Sidechain.h"
#include "Atom.h"
#include "Knotter.h"

Monomer::Monomer()
{
	_backbone = BackbonePtr(new Backbone());
	_sidechain = SidechainPtr(new Sidechain());
}

void Monomer::setup()
{
	_backbone->setMonomer(shared_from_this());
	_sidechain->setMonomer(shared_from_this());
}

void Monomer::addAtom(AtomPtr atom)
{
	_atoms.push_back(atom);
	atom->setMonomer(shared_from_this());
	

	bool isBoth = atom->isBackboneAndSidechain();

	if (isBoth)
	{
		_backbone->addAtom(atom);
		_sidechain->addAtom(atom);
		return;
	}

	bool isBackbone = atom->isBackbone();

	if (isBackbone)
	{
		_backbone->addAtom(atom);
	}
	else
	{
		_sidechain->addAtom(atom);
	}
}


void Monomer::tieAtomsUp()
{
	KnotterPtr knotter = KnotterPtr(new Knotter());
	knotter->setSidechain(_sidechain);

	knotter->tie();
}
