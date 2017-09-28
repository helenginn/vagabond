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
#include "Polymer.h"
#include "Shouter.h"

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

void Monomer::setKick(double value, bool beforeAnchor)
{
	AtomPtr atom;

	if (beforeAnchor)
	{
		atom = getBackbone()->findAtom("C");
	}
	else
	{
		atom = getBackbone()->findAtom("CA");
	}

	ModelPtr model = atom->getModel();
	if (model->getClassName() != "Bond")
	{
		shout_at_helen("Can't kick something that\n isn't a bond!");
	}

	BondPtr bond = ToBondPtr(model);

	Bond::setTorsionBlur(&*bond, value);
}

double Monomer::getKick()
{
	AtomPtr ca = getBackbone()->findAtom("CA");
	ModelPtr model = ca->getModel();
	if (model->getClassName() != "Bond")
	{
		shout_at_helen("Can't kick something that\n isn't a bond!");
	}

	BondPtr bond = ToBondPtr(model);
	
	return Bond::getTorsionBlur(&*bond);
}

void Monomer::tieAtomsUp()
{
	KnotterPtr knotter = KnotterPtr(new Knotter());

	const int start = getPolymer()->getAnchor();

	knotter->setBackbone(_backbone);

	if (getResidueNum() >= start)
	{
		knotter->tieTowardsCTerminus();
	}
	else
	{
		knotter->tieTowardsNTerminus();
	}

	knotter->setSidechain(_sidechain);
	knotter->tie();

	if (getResidueNum() == start)
	{
		BondPtr bond = ToBondPtr(getBackbone()->findAtom("CA")->getModel());
		Bond::setTorsionBlur(&*bond, 0.10);
	}
	else if (getResidueNum() == start - 1)
	{
		BondPtr bond = ToBondPtr(getBackbone()->findAtom("C")->getModel());
		Bond::setTorsionBlur(&*bond, 0.10);
	}

	_backbone->setTied();
	_sidechain->setTied();
}

void Monomer::setConstantDampening(double value)
{
	for (int i = 0; i < modelCount(); i++)
	{
		if (model(i)->getClassName() == "Bond")
		{
			BondPtr bond = ToBondPtr(model(i));
			Bond::setDampening(&*bond, value);
		}
	}
}

bool Monomer::isAfterAnchor()
{
	int anchor = getPolymer()->getAnchor();
	bool after = (_residueNum >= anchor);

	return after;
}
