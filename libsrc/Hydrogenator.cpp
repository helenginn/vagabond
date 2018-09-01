//
// Hydrogenator.cpp
// vagabond
//
// Created by Helen Ginn on 22/03/2018
// Copyright (c) 2018 Helen Ginn
//

#include "Hydrogenator.h"
#include "Element.h"
#include "Bond.h"
#include "Monomer.h"
#include "Backbone.h"
#include "Sidechain.h"
#include "Polymer.h"

Hydrogenator::Hydrogenator()
{
	
}

AtomPtr Hydrogenator::prepareNewHydrogen(AtomPtr parent)
{
	ElementPtr hydrogenElement = Element::getElement("H");

	AtomPtr hydrogen = AtomPtr(new Atom());

	PolymerPtr poly = _monomer->getPolymer();
	
	hydrogen->setFromPDB(false);
	hydrogen->setInitialBFactor(parent->getInitialBFactor());
	hydrogen->setElement(hydrogenElement);
	hydrogen->setOriginalOccupancy(1.);
	hydrogen->setAtomNum(poly->issueAtomNumber());
	_monomer->addAtom(hydrogen);
	
	return hydrogen;
}

bool Hydrogenator::hasHydrogens(BondPtr bond)
{
	for (int i = 0; i < bond->downstreamAtomGroupCount(); i++)
	{
		for (int j = 0; j < bond->downstreamAtomCount(i); j++)
		{
			AtomPtr atom = bond->downstreamAtom(i, j);
			
			if (atom->getElement()->electronCount() == 1)
			{
				return true;
			}
		}
	}
	
	return false;
}

void Hydrogenator::setSpin(AtomList group)
{
	for (int i = 0; i < group.size(); i++)
	{
		AtomPtr atom = group[i].lock();
		
		BondPtr bond = ToBondPtr(atom->getModel());
		BondPtr parent = ToBondPtr(bond->getParentModel());

		bond->setUsingTorsion(true);
		Bond::setKick(&*parent, 2.);
	}
}

void Hydrogenator::setNewGeometry(AtomList group, double bondAngle, 
                                  double torsion, double portion)
{
	for (int i = 0; i < group.size(); i++)
	{
		AtomPtr atom = group[i].lock();
		BondPtr bond = ToBondPtr(atom->getModel());
		
		Bond::setBendAngle(&*bond, deg2rad(bondAngle));	
		BondPtr parent = ToBondPtr(bond->getParentModel());

		if (portion > 0)
		{
			Bond::setCirclePortion(&*bond, portion * 2 * M_PI);
			portion = Bond::getCirclePortion(&*bond);
		}
		
		/* Don't set a torsion angle except for the first atom */
		if (parent->downstreamAtomNum(atom, NULL) > 0)
		{
			continue;
		}
		
		Bond::setTorsion(&*parent, deg2rad(torsion));

	}
}

void Hydrogenator::addHydrogens(AtomList group, int hNum, ...)
{
	va_list arguments;
	va_start(arguments, hNum);
	std::vector<std::string> hNames;
	
	for (int i = 0; i < hNum; i++)
	{
		char *value = va_arg(arguments, char *);
		hNames.push_back(std::string(value));
	}
	
	for (int i = 0; i < group.size(); i++)
	{
		addHydrogens(group[i].lock(), hNames);
	}
}

void Hydrogenator::addHydrogens(AtomPtr minor, std::vector<std::string> hNames)
{
	if (!minor)
	{
		return;
	}
	
	ModelPtr model = minor->getModel();
	
	if (!model) return;
	if (!model->isBond()) return;
	
	BondPtr bond = ToBondPtr(model);
	
	if (hasHydrogens(bond))
	{
		return;
	}
	
	for (int i = 0; i < bond->downstreamAtomGroupCount(); i++)
	{
		/* Find the fraction of the complete "torsion circle" made by the
		* final atom in the downstream atoms. */
		int currentTotal = bond->downstreamAtomCount(i);
		
		int finalTotal = (currentTotal + hNames.size());
		
		double bondAngle = 0;
		
		switch (finalTotal)
		{
			/* Linear */
			case 1:
			bondAngle = deg2rad(180.);
			break;
			
			/* Trigonal */
			case 2:
			bondAngle = deg2rad(120.);
			break;
			
			/* Tetrahedral */
			case 3:
			bondAngle = deg2rad(109.5);
			break;
			
			default:
			break;	
		}
		
		double circlePortion = 0;
		bool hasBonds = bond->downstreamAtomCount(i) > 0;
		
		if (currentTotal > 0)
		{
			circlePortion = bond->getCirclePortion(i, currentTotal - 1);
		}

		if (circlePortion < 0) circlePortion += 1;

		/* Divide the remainder into an appropriate addition per hydrogen. */
		double remaining = 1 - circlePortion;
		
		if (circlePortion > 0.5)
		{
			remaining = -circlePortion;	
		}
		
		int add = currentTotal > 0 ? 1 : 0;
		remaining /= (double)(hNames.size() + add);
		
		if (circlePortion > 0 && false)
		{
			std::cout << bond->description() << std::endl;
			std::cout << "Last circle portion: " << circlePortion << std::endl;
			std::cout << "Adding each time: " << remaining << std::endl;
		}
		
		double nextPortion = circlePortion + remaining;

		for (int j = 0; j < hNames.size(); j++)
		{
			AtomPtr hydrogen = prepareNewHydrogen(minor);
			hydrogen->setAtomName(hNames[j]);

			/* Set the bond length for the new hydrogen */
			BondPtr newBond = BondPtr(new Bond(minor, hydrogen, i));
			newBond->activate();
			Bond::setBondLength(&*newBond, 0.968);

			/* Bond angle... no idea so just going for a tetrahedral thingy */
			Bond::setBendAngle(&*newBond, bondAngle);
			
			/* Set circle portion and increment for the next hydrogen */
			Bond::setCirclePortion(&*newBond, nextPortion * 2 * M_PI);
			nextPortion += remaining;
		}
	}
}

void setHBonds(AtomList list)
{
	for (int i = 0; i < list.size(); i++)
	{
		list[i].lock()->setHBonding(true);
	}
}

void Hydrogenator::hydrogenate()
{
	if (!_monomer)
	{
		return;
	}
	
	BackbonePtr bone = _monomer->getBackbone();

	if (!bone)
	{
		return;	
	}
	
	AtomList nitrogen = bone->findAtoms("N");
	addHydrogens(nitrogen, 1, "H");
	
	/* If anchored to C-terminus, required angle is different */
	AtomPtr prevMajor = bone->betaCarbonTorsionAtom();
	
	double angle = (prevMajor->getAtomName() == "C") ? 124.29 : 117.0;

	setNewGeometry(_monomer->findAtoms("H"), angle, 0.);

	AtomList cAlpha = bone->findAtoms("CA");
	AtomList nitrogens = bone->findAtoms("N");
	AtomList carbonyls = bone->findAtoms("O");
	
	setHBonds(carbonyls);

	if (_monomer->getIdentifier() != "pro")
	{
		setHBonds(nitrogens);
	}
	
	if (_monomer->getIdentifier() == "gly")
	{
		addHydrogens(cAlpha, 2, "HA2", "HA3");
	}
	else
	{
		addHydrogens(cAlpha, 1, "HA");
	}
	
	SidechainPtr side = _monomer->getSidechain();
	
	if (_monomer->getIdentifier() == "met")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		addHydrogens(side->findAtoms("CG"), 2, "HG2", "HG3");
		addHydrogens(side->findAtoms("CE"), 3, "HE1", "HE2", "HE3");
		setNewGeometry(_monomer->findAtoms("HE1"), 109.5, 0., 0.);
		setNewGeometry(_monomer->findAtoms("HE2"), 109.5, 0., 1./3.);
		setNewGeometry(_monomer->findAtoms("HE3"), 109.5, 0., 2./3.);
	}

	if (_monomer->getIdentifier() == "arg")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		addHydrogens(side->findAtoms("CG"), 2, "HG2", "HG3");
		addHydrogens(side->findAtoms("CD"), 2, "HD2", "HD3");
		addHydrogens(side->findAtoms("NE"), 1, "HE");
		addHydrogens(side->findAtoms("NH1"), 2, "HH11", "HH12");
		addHydrogens(side->findAtoms("NH2"), 2, "HH21", "HH22");
		setHBonds(side->findAtoms("NE"));
		setHBonds(side->findAtoms("NH1"));
		setHBonds(side->findAtoms("NH2"));
		setNewGeometry(_monomer->findAtoms("HH11"), 120., 0., 0.);
		setNewGeometry(_monomer->findAtoms("HH12"), 120, 180., 1./2.);
		setNewGeometry(_monomer->findAtoms("HH21"), 120., 0., 0.);
		setNewGeometry(_monomer->findAtoms("HH22"), 120, 180., 1./2.);
	}
	
	if (_monomer->getIdentifier() == "lys")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		addHydrogens(side->findAtoms("CG"), 2, "HG2", "HG3");
		addHydrogens(side->findAtoms("CD"), 2, "HD2", "HD3");
		addHydrogens(side->findAtoms("CE"), 2, "HE1", "HE2");
		addHydrogens(side->findAtoms("NZ"), 3, "HZ1", "HZ2", "HZ3");
		setHBonds(side->findAtoms("NZ"));
	}

	if (_monomer->getIdentifier() == "pro")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		addHydrogens(side->findAtoms("CG"), 2, "HG2", "HG3");
		addHydrogens(side->findAtoms("CD"), 2, "HD2", "HD3");
		setNewGeometry(_monomer->findAtoms("HD2"), 109.5, 140., 0.);
		setNewGeometry(_monomer->findAtoms("HD3"), 109.5, 0., 1./3.);
	}
	
	if (_monomer->getIdentifier() == "ser")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		addHydrogens(side->findAtoms("OG"), 1, "HG");
		setHBonds(side->findAtoms("OG"));
		
		setNewGeometry(_monomer->findAtoms("HG"), 120., 180.);
		setSpin(_monomer->findAtoms("HG"));
	}

	if (_monomer->getIdentifier() == "cys")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		addHydrogens(side->findAtoms("SG"), 1, "HG");

		setNewGeometry(_monomer->findAtoms("HG"), 120., 180.);
		setSpin(_monomer->findAtoms("HG"));
	}
	
	if (_monomer->getIdentifier() == "val")
	{
		addHydrogens(side->findAtoms("CB"), 1, "HB");
		addHydrogens(side->findAtoms("CG1"), 3, "HG11", "HG12", "HG13");
		addHydrogens(side->findAtoms("CG2"), 3, "HG21", "HG22", "HG23");
		setNewGeometry(_monomer->findAtoms("HG11"), 109.5, 0., 0.);
		setNewGeometry(_monomer->findAtoms("HG12"), 109.5, 0., 1./3.);
		setNewGeometry(_monomer->findAtoms("HG13"), 109.5, 0., 2./3.);
		setNewGeometry(_monomer->findAtoms("HG21"), 109.5, 0., 0.);
		setNewGeometry(_monomer->findAtoms("HG22"), 109.5, 0., 1./3.);
		setNewGeometry(_monomer->findAtoms("HG23"), 109.5, 0., 2./3.);
	}
	
	if (_monomer->getIdentifier() == "ala")
	{
		addHydrogens(side->findAtoms("CB"), 3, "HB1", "HB2", "HB3");
		setNewGeometry(_monomer->findAtoms("HB1"), 109.5, 0., 0.);
		setNewGeometry(_monomer->findAtoms("HB2"), 109.5, 0., 1./3.);
		setNewGeometry(_monomer->findAtoms("HB3"), 109.5, 0., 2./3.);
	}

	if (_monomer->getIdentifier() == "his")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		addHydrogens(side->findAtoms("ND1"), 1, "HD1");
		addHydrogens(side->findAtoms("CE1"), 1, "HE1");
		addHydrogens(side->findAtoms("CD2"), 1, "HD2");
		addHydrogens(side->findAtoms("NE2"), 1, "HE2");
		
		setNewGeometry(_monomer->findAtoms("HE1"), 120., 180.);
		setNewGeometry(_monomer->findAtoms("HE2"), 120., 180.);
	}

	if (_monomer->getIdentifier() == "tyr")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		addHydrogens(side->findAtoms("CD1"), 1, "HD1");
		addHydrogens(side->findAtoms("CE1"), 1, "HE1");
		addHydrogens(side->findAtoms("CE2"), 1, "HE2");
		addHydrogens(side->findAtoms("CD2"), 1, "HD2");
		addHydrogens(side->findAtoms("OH"), 1, "HH");
		setHBonds(side->findAtoms("OH"));
		
		setNewGeometry(_monomer->findAtoms("HE1"), 120., 180.);
		setNewGeometry(_monomer->findAtoms("HH"), 120., 180.);
		setSpin(_monomer->findAtoms("HH"));
	}

	if (_monomer->getIdentifier() == "phe")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		addHydrogens(side->findAtoms("CD1"), 1, "HD1");
		addHydrogens(side->findAtoms("CE1"), 1, "HE1");
		addHydrogens(side->findAtoms("CE2"), 1, "HE2");
		addHydrogens(side->findAtoms("CD2"), 1, "HD2");
		addHydrogens(side->findAtoms("CZ"), 1, "HZ");
		
		setNewGeometry(_monomer->findAtoms("HE1"), 120., 180.);
		setNewGeometry(_monomer->findAtoms("HZ"), 120., 180.);

	}

	if (_monomer->getIdentifier() == "trp")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		addHydrogens(side->findAtoms("CD1"), 1, "HD1");
		addHydrogens(side->findAtoms("NE1"), 1, "HE1");
		addHydrogens(side->findAtoms("CE3"), 1, "HE3");
		addHydrogens(side->findAtoms("CZ2"), 1, "HZ2");
		addHydrogens(side->findAtoms("CZ3"), 1, "HZ3");
		addHydrogens(side->findAtoms("CH2"), 1, "HH2");
		
		setNewGeometry(_monomer->findAtoms("HH2"), 120., 180.);
		setNewGeometry(_monomer->findAtoms("HZ2"), 120., 0.);

	}
	
	if (_monomer->getIdentifier() == "ile")
	{
		addHydrogens(side->findAtoms("CB"), 1, "HB");
		addHydrogens(side->findAtoms("CG1"), 2, "HG12", "HG13");
		addHydrogens(side->findAtoms("CD1"), 3, "HD11", "HD12", "HD13");
		addHydrogens(side->findAtoms("CG2"), 3, "HG21", "HG22", "HG23");
		setNewGeometry(_monomer->findAtoms("HG21"), 109.5, 0., 0.);
		setNewGeometry(_monomer->findAtoms("HG22"), 109.5, 0., 1./3.);
		setNewGeometry(_monomer->findAtoms("HG23"), 109.5, 0., 2./3.);
		setNewGeometry(_monomer->findAtoms("HD11"), 109.5, 0., 0.);
		setNewGeometry(_monomer->findAtoms("HD12"), 109.5, 0., 1./3.);
		setNewGeometry(_monomer->findAtoms("HD13"), 109.5, 0., 2./3.);
	}
	
	if (_monomer->getIdentifier() == "leu")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		addHydrogens(side->findAtoms("CG"), 1, "HG");
		addHydrogens(side->findAtoms("CD1"), 3, "HD11", "HD12", "HD13");
		addHydrogens(side->findAtoms("CD2"), 3, "HD21", "HD22", "HD23");
		setNewGeometry(_monomer->findAtoms("HD11"), 109.5, 0., 0.);
		setNewGeometry(_monomer->findAtoms("HD12"), 109.5, 0., 1./3.);
		setNewGeometry(_monomer->findAtoms("HD13"), 109.5, 0., 2./3.);
		setNewGeometry(_monomer->findAtoms("HD21"), 109.5, 0., 0.);
		setNewGeometry(_monomer->findAtoms("HD22"), 109.5, 0., 1./3.);
		setNewGeometry(_monomer->findAtoms("HD23"), 109.5, 0., 2./3.);
	}

	if (_monomer->getIdentifier() == "asp")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		setHBonds(side->findAtoms("OD1"));
		setHBonds(side->findAtoms("OD2"));
	}
	
	if (_monomer->getIdentifier() == "asn")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		addHydrogens(side->findAtoms("ND2"), 2, "HD21", "HD22");
		setHBonds(side->findAtoms("OD1"));
		setHBonds(side->findAtoms("ND2"));
		setNewGeometry(_monomer->findAtoms("HD21"), 120., 0., 0.);
		setNewGeometry(_monomer->findAtoms("HD22"), 120, 180., 1./2.);
	}

	if (_monomer->getIdentifier() == "glu")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		addHydrogens(side->findAtoms("CG"), 2, "HG2", "HG3");
		setHBonds(side->findAtoms("OD1"));
		setHBonds(side->findAtoms("OD2"));
	}

	if (_monomer->getIdentifier() == "gln")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		addHydrogens(side->findAtoms("CG"), 2, "HG2", "HG3");
		addHydrogens(side->findAtoms("NE2"), 2, "HE21", "HE22");
		setHBonds(side->findAtoms("OD1"));
		setHBonds(side->findAtoms("ND2"));

		setNewGeometry(_monomer->findAtoms("HE21"), 120., 0., 0.);
		setNewGeometry(_monomer->findAtoms("HE22"), 120, 180., 1./2.);
	}

	if (_monomer->getIdentifier() == "thr")
	{
		addHydrogens(side->findAtoms("CB"), 1, "HB");
		addHydrogens(side->findAtoms("OG1"), 1, "HG1");
		addHydrogens(side->findAtoms("CG2"), 3, "HG21", "HG22", "HG23");
		setHBonds(side->findAtoms("OG1"));

		setNewGeometry(_monomer->findAtoms("HG1"), 120, 180., 0.);
		setNewGeometry(_monomer->findAtoms("HG21"), 109.5, 0., 0.);
		setNewGeometry(_monomer->findAtoms("HG22"), 109.5, 0., 1./3.);
		setNewGeometry(_monomer->findAtoms("HG23"), 109.5, 0., 2./3.);
	}
}














