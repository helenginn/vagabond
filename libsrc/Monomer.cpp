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
#include "Bond.h"
#include "shared_ptrs.h"

Monomer::Monomer()
{
	_backbone = BackbonePtr(new Backbone());
	_sidechain = SidechainPtr(new Sidechain());
}

void Monomer::setup()
{
	_backbone->setMonomer(shared_from_this());
	_backbone->setPolymer(getPolymer());
	_backbone->setResNum(_residueNum);
	_sidechain->setMonomer(shared_from_this());
	_sidechain->setResNum(_residueNum);

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

	const int start = 85;

	if (getResidueNum() >= start && getResidueNum() <= 124)
	{
		bool useAbsolute = (getResidueNum() < start + 12);

		if (useAbsolute)
		{
			getBackbone()->setUseAbsolute();
			getSidechain()->setUseAbsolute();
		}
		knotter->setBackbone(_backbone);
		knotter->tieTowardsCTerminus();
		knotter->setSidechain(_sidechain);
		knotter->tie();

		if (getResidueNum() == start)
		{
			BondPtr bond = ToBondPtr(getBackbone()->findAtom("CA")->getModel());
			Bond::setTorsionBlur(&*bond, 0.0);
		}

		_backbone->setTied();
		_sidechain->setTied();
	}
}
