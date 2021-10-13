//
// Hydrogenator.cpp
// vagabond
//
// Created by Helen Ginn on 22/03/2018
// Copyright (c) 2018 Helen Ginn
//

#include "Hydrogenator.h"
#include "Bond.h"
#include "Element.h"
#include "Monomer.h"
#include "Backbone.h"
#include "Sidechain.h"
#include "Polymer.h"
#include "Options.h"
#include "Crystal.h"

Hydrogenator::Hydrogenator()
{
	_cAlpha2H = 0.94;
	
}

AtomPtr Hydrogenator::prepareNewHydrogen(AtomPtr parent)
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	ElementPtr hydrogenElement = Element::getElement("H");

	AtomPtr hydrogen = AtomPtr(new Atom());

	PolymerPtr poly = _monomer->getPolymer();
	
	hydrogen->setFromPDB(false);
	hydrogen->setInitialBFactor(parent->getInitialBFactor());
	hydrogen->setElement(hydrogenElement);
	hydrogen->setOriginalOccupancy(1.);
	hydrogen->setAtomNum(crystal->issueAtomNumber());
	
	return hydrogen;
}

bool Hydrogenator::hasHydrogens(BondPtr bond)
{
	for (int i = 0; i < bond->downstreamBondGroupCount(); i++)
	{
		for (int j = 0; j < bond->downstreamBondCount(i); j++)
		{
			AtomPtr atom = bond->downstreamAtom(i, j);
			
			if (atom->getElectronCount() == 1)
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
		AtomPtr atom = group[i];
		
		if (atom->getModel() && atom->getModel()->isBond())
		{
			BondPtr bond = ToBondPtr(atom->getModel());
			bond->setUsingTorsion(true);
		}
		
	}
}

void Hydrogenator::setNewGeometry(AtomList group, double bondAngle, 
                                  double torsion, double portion)
{
	for (int i = 0; i < group.size(); i++)
	{
		AtomPtr atom = group[i];
		
		if (!atom->getModel()->isBond())
		{
			continue;
		}
		
		BondPtr bond = ToBondPtr(atom->getModel());
		
		Bond::setBendAngle(&*bond, deg2rad(bondAngle));	

		if (portion > 0)
		{
			Bond::setCirclePortion(&*bond, portion * 2 * M_PI);
			portion = Bond::getCirclePortion(&*bond);
		}

		adjustBond(bond);
		
		ModelPtr mod = bond->getParentModel();
		
		if (!mod->isBond())
		{
			continue;
		}
		
		BondPtr parent = ToBondPtr(mod);
		int num = parent->downstreamBondNum(&*bond, NULL);
		
		/* Don't set a torsion angle except for the first atom */
		if (num > 0 || torsion < 0)
		{
			continue;
		}
		
		Bond::setTorsion(&*bond, deg2rad(torsion));
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
		addHydrogens(group[i], hNames);
	}
}

void Hydrogenator::adjustProlineHydrogens(BondPtr newBond, bool current)
{
	AtomPtr minor = newBond->getMajor();
	BondPtr nBond = ToBondPtr(minor->getModel());
	AtomPtr major = nBond->getMajor();

	if (((newBond->getMinor()->getAtomName() == "HD2") ||
	     (newBond->getMinor()->getAtomName() == "HD3")) &&
		(minor->getAtomName() == "CD") &&
		(minor->getMonomer()->getResCode() == "P"))
	{
		AtomPtr n = minor->getMonomer()->findAtom("N");
		AtomPtr cg = major;
		AtomPtr cb = ToBondPtr(major->getModel())->getMajor();
		AtomPtr cd = minor;
		vec3 pos_n = n->getInitialPosition();
		vec3 pos_cd = cd->getInitialPosition();
		vec3 pos_cg = cg->getInitialPosition();
		vec3 pos_cb = cb->getInitialPosition();
		
		if (current)
		{
			pos_n = n->getAbsolutePosition();
			pos_cd = cd->getAbsolutePosition();
			pos_cg = cg->getAbsolutePosition();
			pos_cb = cb->getAbsolutePosition();
		}

		double ratio = newBond->getGeomRatio();
		double angleh2, angleh3;
		mat3x3 m = ExplicitModel::makeTorsionBasis(pos_cg, pos_n, 
		                                           pos_cd, empty_vec3());
		vec3 pos_h2 = Bond::positionFromTorsion(m, deg2rad(-120),
		                                       0.968, ratio, pos_cd);
		ExplicitModel::makeTorsionBasis(pos_cb, pos_cg, pos_cd, pos_h2, 
		                                &angleh2);

		vec3 pos_h3 = Bond::positionFromTorsion(m, deg2rad(120),
		                                       0.968, ratio, pos_cd);
		ExplicitModel::makeTorsionBasis(pos_cb, pos_cg, pos_cd, pos_h3, 
		                                &angleh3);

		if (newBond->getMinor()->getAtomName() == "HD2")
		{
			Bond::setTorsion(&*newBond, angleh2);
		}
		if (newBond->getMinor()->getAtomName() == "HD3")
		{
			Bond::setCirclePortion(&*newBond, angleh3 - angleh2);
		}
	}
}

void Hydrogenator::adjustBond(BondPtr newBond)
{
	AtomPtr minor = newBond->getMajor();
	BondPtr nBond = ToBondPtr(minor->getModel());
	AtomPtr major = nBond->getMajor();

	if ((minor->getAtomName() == "N") &&
	    (major->getAtomName() == "CA"))
	{
		Bond::setBendAngle(&*newBond, deg2rad(114.0));
	}
	else if ((minor->getAtomName() == "N") &&
	    (major->getAtomName() == "C"))
	{
		Bond::setBendAngle(&*newBond, deg2rad(123.9));
	}
	else if ((minor->getAtomName() == "CA") &&
	         (major->getAtomName() == "C"))
	{
		Bond::setBendAngle(&*newBond, deg2rad(108.0));
	}

	adjustProlineHydrogens(newBond, false);
}

double Hydrogenator::getHBondLength(AtomPtr minor)
{
	std::string id = minor->getMonomer()->getIdentifier();
	std::string name = minor->getAtomName();

	if ((id == "tyr" || id == "phe" || id == "his") && 
	    (name == "CD1" || name == "CD2" || name == "CE1" || name == "CE2"
	    || name == "CZ"))
	{
		return 0.94;
	}

	if (minor->getElement()->getSymbol() == "N")
	{
		if (minor->getAtomName() == "N")
		{
//			return 0.911;
		}
		return 0.86;
	}
	else if (minor->getElement()->getSymbol() == "O")
	{
		return 0.846;
	}
	else if (name == "CA")
	{
		return _cAlpha2H;
	}
	
	return 0.968;
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

	/* Find the fraction of the complete "torsion circle" made by the
	 * final atom in the downstream atoms. */
	int currentTotal = 0;
	
	if (bond->downstreamBondGroupCount() > 0)
	{
		currentTotal = bond->downstreamBondCount(0);
	}

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
	std::string alt = "";
	alt = bond->getMinor()->getAlternativeConformer();

	/* We want to add onto the existing circle portion */
	if (currentTotal > 0)
	{
		/* Returned in radians */
		double portion = 0;
		BondPtr last = bond->downstreamBond(0, currentTotal - 1);
		circlePortion = Bond::getCirclePortion(&*last);
	}

	if (circlePortion < 0) circlePortion += deg2rad(360);

	/* Divide the remainder into an appropriate addition per hydrogen. */
	double remaining = deg2rad(360) - circlePortion;

	if (circlePortion > 0.5 * deg2rad(360))
	{
		remaining = -circlePortion;	
	}

	/* We don't need to add to the first hydrogen if we don't
	 * have any original non-hydrogen downstream bonds */
	int add = currentTotal > 0 ? 1 : 0;
	remaining /= (double)(hNames.size() + add);

	double nextPortion = circlePortion + remaining * add;
	CrystalPtr crystal = Options::getActiveCrystal();

	for (int j = 0; j < hNames.size(); j++)
	{
		AtomPtr hydrogen = prepareNewHydrogen(minor);
		hydrogen->setAtomName(hNames[j]);
		_monomer->addAtom(hydrogen);
		crystal->addAtom(hydrogen);
		hydrogen->setAlternativeConformer(alt);

		/* Set the bond length for the new hydrogen */
		BondPtr newBond = BondPtr(new Bond(minor, hydrogen, 0));
		newBond->activate();
		double length = getHBondLength(minor);
		Bond::setBondLength(&*newBond, length);

		/* Bond angle... no idea so just going for a tetrahedral thingy */
		Bond::setBendAngle(&*newBond, bondAngle);

		/* Set circle portion and increment for the next hydrogen */
		Bond::setCirclePortion(&*newBond, nextPortion);
		nextPortion += remaining;
	}
}

void setHBonds(AtomList list)
{
	for (int i = 0; i < list.size(); i++)
	{
		list[i]->setHBonding(true);
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
	if (_monomer->getIdentifier() != "pro")
	{
		addHydrogens(nitrogen, 1, "H");
	}
	
	/* If anchored to C-terminus, required angle is different */
	AtomPtr prevMajor = bone->betaCarbonTorsionAtom();
	if (!prevMajor)
	{
		return;
	}
	
	double angle = (prevMajor->getAtomName() == "C") ? 124.29 : 117.0;

	setNewGeometry(_monomer->findAtoms("H"), angle);

	AtomList cAlpha = bone->findAtoms("CA");
	AtomList nitrogens = bone->findAtoms("N");
	AtomList carbonyls = bone->findAtoms("O");
	
	setHBonds(carbonyls);

	setHBonds(nitrogens);
	
	if (_monomer->getIdentifier() == "gly")
	{
		addHydrogens(cAlpha, 2, "HA2", "HA3");
		setNewGeometry(_monomer->findAtoms("HA2"), 109.5);
		setNewGeometry(_monomer->findAtoms("HA3"), 109.5);
	}
	else
	{
		addHydrogens(cAlpha, 1, "HA");
		setNewGeometry(_monomer->findAtoms("HA"), 109.5);
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
		setNewGeometry(_monomer->findAtoms("HD2"), 112.5, -1, 0.);
		setNewGeometry(_monomer->findAtoms("HD3"), 112.5, 0., 1./3.);
	}
	
	if (_monomer->getIdentifier() == "ser")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		addHydrogens(side->findAtoms("OG"), 1, "HG");
		setHBonds(side->findAtoms("OG"));
		
		setNewGeometry(_monomer->findAtoms("HG"), 109.5, 180.);
		setSpin(_monomer->findAtoms("HG"));
	}

	if (_monomer->getIdentifier() == "cys")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
//		addHydrogens(side->findAtoms("SG"), 1, "HG");

		setSpin(_monomer->findAtoms("HG"));
	}
	
	if (_monomer->getIdentifier() == "val")
	{
		addHydrogens(side->findAtoms("CB"), 1, "HB");
		addHydrogens(side->findAtoms("CG1"), 3, "HG11", "HG12", "HG13");
		addHydrogens(side->findAtoms("CG2"), 3, "HG21", "HG22", "HG23");
		setNewGeometry(_monomer->findAtoms("HG11"), 109.5, 60., 0.);
		setNewGeometry(_monomer->findAtoms("HG12"), 109.5, 0., 1./3.);
		setNewGeometry(_monomer->findAtoms("HG13"), 109.5, 0., 2./3.);
		setNewGeometry(_monomer->findAtoms("HG21"), 109.5, 60., 0.);
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
		
		setNewGeometry(_monomer->findAtoms("HD1"), 126.0, 180.);
		setNewGeometry(_monomer->findAtoms("HD2"), 125.4, 180.);
		setNewGeometry(_monomer->findAtoms("HE1"), 124.8, 180.);
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
		setNewGeometry(_monomer->findAtoms("HG21"), 109.5, 60., 0.);
		setNewGeometry(_monomer->findAtoms("HG22"), 109.5, 0., 1./3.);
		setNewGeometry(_monomer->findAtoms("HG23"), 109.5, 0., 2./3.);
		setNewGeometry(_monomer->findAtoms("HD11"), 109.5, 60., 0.);
		setNewGeometry(_monomer->findAtoms("HD12"), 109.5, 0., 1./3.);
		setNewGeometry(_monomer->findAtoms("HD13"), 109.5, 0., 2./3.);
	}
	
	if (_monomer->getIdentifier() == "leu")
	{
		addHydrogens(side->findAtoms("CB"), 2, "HB2", "HB3");
		addHydrogens(side->findAtoms("CG"), 1, "HG");
		addHydrogens(side->findAtoms("CD1"), 3, "HD11", "HD12", "HD13");
		addHydrogens(side->findAtoms("CD2"), 3, "HD21", "HD22", "HD23");
		setNewGeometry(_monomer->findAtoms("HD11"), 109.5, 60., 0.);
		setNewGeometry(_monomer->findAtoms("HD12"), 109.5, 0., 1./3.);
		setNewGeometry(_monomer->findAtoms("HD13"), 109.5, 0., 2./3.);
		setNewGeometry(_monomer->findAtoms("HD21"), 109.5, 60., 0.);
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

		setNewGeometry(_monomer->findAtoms("HG1"), 109.5, 180., 0.);
		setNewGeometry(_monomer->findAtoms("HG21"), 109.5, 0., 0.);
		setNewGeometry(_monomer->findAtoms("HG22"), 109.5, 0., 1./3.);
		setNewGeometry(_monomer->findAtoms("HG23"), 109.5, 0., 2./3.);
	}
}

