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
#include "Options.h"

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

void Monomer::removeAtom(AtomPtr atom)
{
	_sidechain->removeAtom(atom);
	_backbone->removeAtom(atom);
	
	AtomGroup::removeAtom(atom);
}

void Monomer::empty()
{
	_backbone->empty();
	_sidechain->empty();
	AtomGroup::empty();
}

void Monomer::addAtom(AtomPtr atom)
{
	AtomGroup::addAtom(atom);

	atom->setMonomer(shared_from_this());

	if (getPolymer())
	{
		getPolymer()->addAtom(atom);
	}

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

	const int start = getPolymer()->getAnchor();

	knotter->setBackbone(_backbone);

	try
	{
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
	}
	catch (AtomPtr a)
	{
		printList();

		std::cout << std::endl <<"Backbone: " << std::endl;
		_backbone->printList();

		std::cout << std::endl <<"Sidechain: " << std::endl;
		_sidechain->printList();
		shout_at_user("Could not tie protein due to a problem around\n"\
		              "atom " + a->shortDesc() + "! Please check all "\
		              "your atoms are\nin the model for this monomer.");
		std::cout << "Monomer: " << std::endl;

		exit(1);
	}

	if (getResidueNum() == start)
	{
		BondPtr bond = ToBondPtr(getBackbone()->findAtom("CA")->getModel());
	}
	else if (getResidueNum() == start - 1)
	{
		BondPtr bond = ToBondPtr(getBackbone()->findAtom("C")->getModel());
	}

	_backbone->setTied();
	_sidechain->setTied();
}

void Monomer::refine(CrystalPtr target,
                     RefinementType rType)
{
	BackbonePtr backbone = getBackbone();

	if (backbone)
	{
		backbone->refine(target, rType);
	}


	SidechainPtr victim = getSidechain();

	if (victim && victim->canRefine())
	{
		victim->refine(target, rType);
	}
}

bool Monomer::isAfterAnchor()
{
	int anchor = getPolymer()->getAnchor();
	bool after = (_residueNum >= anchor);

	return after;
}

std::string Monomer::getResCode()
{
	std::string id = getIdentifier();
	return GeomTable::getResCode(id);
}

void Monomer::addProperties()
{
	addStringProperty("identifier", &_identifier);
	addIntProperty("res_num", &_residueNum);
	addChild("sidechain", _sidechain);
	addChild("backbone", _backbone);
	
	AtomGroup::addProperties();
}

void Monomer::addObject(ParserPtr object, std::string category)
{
	if (category == "sidechain")
	{
		SidechainPtr side = ToSidechainPtr(object);
		_sidechain = side; 
		_sidechain->setMonomer(shared_from_this());
	}

	if (category == "backbone")
	{
		BackbonePtr back = ToBackbonePtr(object);
		_backbone = back;
		_backbone->setMonomer(shared_from_this());
	}

	AtomGroup::addObject(object, category);
}

void Monomer::linkReference(BaseParserPtr object, std::string category)
{
	if (category == "atom")
	{
		AtomPtr atom = ToAtomPtr(object);
		addAtom(atom);
	}   
}

void Monomer::postParseTidy()
{

}

void Monomer::setPolymer(PolymerPtr pol)
{
	_myPolymer = pol;
	
	if (_backbone)
	{
		_backbone->setPolymer(pol);
	}
	
	if (_sidechain)
	{
		_sidechain->setPolymer(pol);
	}
}
