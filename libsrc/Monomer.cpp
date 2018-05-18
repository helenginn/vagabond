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

void Monomer::addAtom(AtomPtr atom)
{
	std::string newAtomName = atom->getAtomName();

	AtomPtr existingAtom = findAtom(newAtomName);

	if (existingAtom)
	{
		//    std::string atomDesc = existingAtom->shortDesc();
		//    warn_user("There is an existing conformation for " + atomDesc
		//              + "\nCurrently only side chains are supported.");
	}

	AtomGroup::addAtom(atom);

	atom->setMonomer(shared_from_this());

	if (getPolymer())
	{
		getPolymer()->addAtom(atom);
	}

	bool isBoth = atom->isBackboneAndSidechain();

	if (isBoth)
	{
		if (!existingAtom)
		{
			_backbone->addAtom(atom);
		}

		_sidechain->addAtom(atom);
		return;
	}

	bool isBackbone = atom->isBackbone();

	if (isBackbone && !existingAtom)
	{
		_backbone->addAtom(atom);
	}
	else if (!isBackbone)
	{
		_sidechain->addAtom(atom);
	}
}

void Monomer::setSideKick(double value)
{
	AtomPtr atom = findAtom("CB");

	if (!atom)
	{
		return;
	}

	ModelPtr model = atom->getModel();
	if (!model->isBond())
	{
		return;
	}

	BondPtr bond = ToBondPtr(model);

	if (bond->isFixed() || !bond->isUsingTorsion())
	{
		return;
	}

	Bond::setTorsionBlur(&*bond, value);
}

void Monomer::setKick(double value, bool beforeAnchor)
{
	AtomPtr atom;

	if (beforeAnchor)
	{
		atom = findAtom("C");
	}
	else
	{
		atom = findAtom("CA");
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
	AtomPtr ca = findAtom("CA");
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
	}
	else if (getResidueNum() == start - 1)
	{
		BondPtr bond = ToBondPtr(getBackbone()->findAtom("C")->getModel());
	}

	_backbone->setTied();
	_sidechain->setTied();
}

double Monomer::vsRefine(void *object)
{
	OptionsPtr options = Options::getRuntimeOptions();
	CrystalPtr active = options->getActiveCrystal();

	Parser *parser = static_cast<Parser *>(object);
	Monomer *monomer = dynamic_cast<Monomer *>(parser);

	SidechainPtr victim = monomer->getSidechain();
	
	if (!victim || !victim->canRefine())
	{
		return 0;
	}

	std::cout << "Refining monomer " << monomer->getResidueNum() << std::endl;
	victim->addParamType(ParamOptionTorsion, 0.04);
	victim->addParamType(ParamOptionBondAngle, 0.5);
	victim->addParamType(ParamOptionKick, 0.25);
	victim->addParamType(ParamOptionDampen, 0.10);
	victim->addParamType(ParamOptionMagicAngles, 20.);
	victim->addParamType(ParamOptionNumBonds, 5);

	victim->refine(active, RefinementFine);

	return 0;
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

void Monomer::setBackboneDampening(double value)
{
	for (int i = 0; i < getBackbone()->atomCount(); i++)
	{
		if (getBackbone()->atom(i)->getModel()->isBond())
		{
			ModelPtr model = atom(i)->getModel();
			
			if (!model || !model->isBond()) continue;

			BondPtr bond = ToBondPtr(model);
			
			if (!bond->connectsAtom("CA"))
			{
				continue;
			}

			if (bond->isRefinable())
			{
				Bond::setDampening(&*bond, value);
			}
		}
	}
}

void Monomer::setSidechainDampening(double value)
{
	for (int i = 0; i < getSidechain()->atomCount(); i++)
	{
		ModelPtr model = getSidechain()->atom(i)->getModel();

		if (!model || !model->isBond())
		{
			continue;
		}

		if (ToBondPtr(model)->isRefinable())
		{
			Bond::setDampening(&*model, value);
		}
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
	
	exposeFunction("refine_sidechain_to_density", vsRefine);
	
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

void Monomer::linkReference(ParserPtr object, std::string category)
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
