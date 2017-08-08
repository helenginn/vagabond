//
//  Knotter.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Knotter.h"
#include "Shouter.h"
#include "Absolute.h"
#include "Sidechain.h"
#include "BackBone.h"
#include "Monomer.h"
#include <iostream>
#include "Atom.h"
#include "Bond.h"

void Knotter::tie()
{
	if (!_sidechain)
	{
		shout_at_helen("Helen set up a Knotter to tie up atoms,\n"\
					   "and then proceeded to try to convert them\n"\
					   "without specifying a sidechain first.\n"\
					   "What an idiot. Haha!");
	}

	MonomerPtr monomer = _sidechain->getMonomer();
	std::string residue = monomer->getIdentifier();
	int resNum = monomer->getResidueNum();

	bool convertable = false;

	if (resNum != 86)
	{
		return;
	}

	if (residue == "lys")
	{
		convertable = true;
		makeLysine();
	}

	if (residue == "ser")
	{
		convertable = true;
		makeSerine();
	}

	if (residue == "cys")
	{
		convertable = true;
		makeCysteine();
	}

	if (residue == "val")
	{
		convertable = true;
		makeValine();
	}

	if (residue == "thr")
	{
		convertable = true;
		makeThreonine();
	}

	if (convertable)
	{
		std::cout << "Knotter is tying up residue: " << residue
		<< resNum << std::endl;
	}
	else
	{
//		std::cout << "Knotter doesn't know how to tie up " << residue << std::endl;
	}
}

void Knotter::makeLysine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr nSpine = backbone->findAtom("N");
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr hBackbone = _sidechain->findAtom("HA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr hBeta2 = _sidechain->findAtom("HB2");
	AtomPtr hBeta3 = _sidechain->findAtom("HB3");
	AtomPtr cGamma = _sidechain->findAtom("CG");
	AtomPtr hGamma2 = _sidechain->findAtom("HG2");
	AtomPtr hGamma3 = _sidechain->findAtom("HG3");
	AtomPtr cDelta = _sidechain->findAtom("CD");
	AtomPtr hDelta2 = _sidechain->findAtom("HD2");
	AtomPtr hDelta3 = _sidechain->findAtom("HD3");
	AtomPtr cEpsilon = _sidechain->findAtom("CE");
	AtomPtr hEpsilon2 = _sidechain->findAtom("HE2");
	AtomPtr hEpsilon3 = _sidechain->findAtom("HE3");
	AtomPtr nOmega = _sidechain->findAtom("NZ");
	AtomPtr hOmega1 = _sidechain->findAtom("HZ1");
	AtomPtr hOmega2 = _sidechain->findAtom("HZ2");
	AtomPtr hOmega3 = _sidechain->findAtom("HZ3");

	AbsolutePtr inherit = std::static_pointer_cast<Absolute>(cAlpha->getModel());

	BondPtr ca2cb = BondPtr(new Bond(cAlpha, cBeta));
	ca2cb->setBendTowards(hBackbone);
	ca2cb->setTorsionAtoms(nSpine, cGamma);
	ca2cb->activate(_sidechain, inherit);

	BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
	cb2cg->setTorsionAtoms(cAlpha, cDelta);
	cb2cg->activate(_sidechain, inherit);

	BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
	BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));

	BondPtr cg2cd = BondPtr(new Bond(cGamma, cDelta));
	cg2cd->setTorsionAtoms(cBeta, cEpsilon);
	cg2cd->activate(_sidechain, inherit);

	BondPtr cg2hg2 = BondPtr(new Bond(cGamma, hGamma2));
	BondPtr cg2hg3 = BondPtr(new Bond(cGamma, hGamma3));

	BondPtr cd2ce = BondPtr(new Bond(cDelta, cEpsilon));
	cd2ce->setTorsionAtoms(cGamma, nOmega);
	cd2ce->activate(_sidechain, inherit);

	BondPtr cd2hd2 = BondPtr(new Bond(cDelta, hDelta2));
	BondPtr cd2hd3 = BondPtr(new Bond(cDelta, hDelta3));

	BondPtr ce2nz = BondPtr(new Bond(cEpsilon, nOmega));
	ce2nz->activate(_sidechain, inherit);

	BondPtr ce2he2 = BondPtr(new Bond(cEpsilon, hEpsilon2));
	BondPtr ce2he3 = BondPtr(new Bond(cEpsilon, hEpsilon3));

	BondPtr nz2hz1 = BondPtr(new Bond(nOmega, hOmega1));
	BondPtr nz2hz2 = BondPtr(new Bond(nOmega, hOmega2));
	BondPtr nz2hz3 = BondPtr(new Bond(nOmega, hOmega3));

	nz2hz1->activate(_sidechain, inherit);
	nz2hz2->activate(_sidechain, inherit);
	nz2hz3->activate(_sidechain, inherit);

	ce2he2->activate(_sidechain, inherit);
	ce2he3->activate(_sidechain, inherit);

	cd2hd2->activate(_sidechain, inherit);
	cd2hd3->activate(_sidechain, inherit);

	cg2hg2->activate(_sidechain, inherit);
	cg2hg3->activate(_sidechain, inherit);

	cb2hb3->activate(_sidechain, inherit);
	cb2hb2->activate(_sidechain, inherit);
}

void Knotter::makeSerine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);
	//int resNum = monomer->getResidueNum();

	AtomPtr nSpine = backbone->findAtom("N");
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr hBeta2 = _sidechain->findAtom("HB2");
	AtomPtr hBeta3 = _sidechain->findAtom("HB3");
	AtomPtr oGamma = _sidechain->findAtom("OG");
	AtomPtr hGamma = _sidechain->findAtom("HG");
	AtomPtr hBackBone = _sidechain->findAtom("HA");

	AbsolutePtr inherit = std::static_pointer_cast<Absolute>(cAlpha->getModel());

	if (!hBeta2 || !hBeta3 || !hGamma)
	{
		warn_user("Missing some hydrogens in " + residue + ", continuing.");
	}

	if (!cAlpha || !cBeta || !oGamma)
	{
		warn_user("Missing vital atoms in " + residue + ", giving up.");
		return;
	}

	BondPtr ca2cb = BondPtr(new Bond(cAlpha, cBeta));
//	Bond::setBondLength(&*ca2cb, 1.530);
	ca2cb->setBendTowards(hBackBone);
	ca2cb->setTorsionAtoms(nSpine, oGamma);
	ca2cb->activate(_sidechain, inherit);

	BondPtr cb2og = BondPtr(new Bond(cBeta, oGamma));
	cb2og->activate(_sidechain, inherit);

	BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
	cb2hb2->activate(_sidechain, inherit);
	BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));
	cb2hb3->activate(_sidechain, inherit);

	BondPtr og2hg = BondPtr(new Bond(oGamma, hGamma));
	og2hg->activate(_sidechain, inherit);
}

void Knotter::makeCysteine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr nSpine = backbone->findAtom("N");
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr hBeta2 = _sidechain->findAtom("HB2");
	AtomPtr hBeta3 = _sidechain->findAtom("HB3");
	AtomPtr sGamma = _sidechain->findAtom("SG");
	AtomPtr hGamma = _sidechain->findAtom("HG");
	AtomPtr hBackBone = _sidechain->findAtom("HA");

	AbsolutePtr inherit = std::static_pointer_cast<Absolute>(cAlpha->getModel());

	BondPtr ca2cb = BondPtr(new Bond(cAlpha, cBeta));
	Bond::setBondLength(&*ca2cb, 1.530);
	ca2cb->setBendTowards(hBackBone);
	ca2cb->setTorsionAtoms(nSpine, sGamma);
	ca2cb->activate(_sidechain, inherit);

	BondPtr cb2sg = BondPtr(new Bond(cBeta, sGamma));
	cb2sg->activate(_sidechain, inherit);

	BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
	cb2hb2->activate(_sidechain, inherit);
	BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));
	cb2hb3->activate(_sidechain, inherit);

	BondPtr sg2hg = BondPtr(new Bond(sGamma, hGamma));
	sg2hg->activate(_sidechain, inherit);
}

void Knotter::makeValine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr nSpine = backbone->findAtom("N");
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr hBeta = _sidechain->findAtom("HB");
	AtomPtr cGamma1 = _sidechain->findAtom("CG1");
	AtomPtr hGamma11 = _sidechain->findAtom("HG11");
	AtomPtr hGamma12 = _sidechain->findAtom("HG12");
	AtomPtr hGamma13 = _sidechain->findAtom("HG13");
	AtomPtr cGamma2 = _sidechain->findAtom("CG2");
	AtomPtr hGamma21 = _sidechain->findAtom("HG21");
	AtomPtr hGamma22 = _sidechain->findAtom("HG22");
	AtomPtr hGamma23 = _sidechain->findAtom("HG23");
	AtomPtr hBackbone = _sidechain->findAtom("HA");

	AbsolutePtr inherit = std::static_pointer_cast<Absolute>(cAlpha->getModel());

	BondPtr ca2cb = BondPtr(new Bond(cAlpha, cBeta));
	ca2cb->setTorsionAtoms(nSpine, cGamma1);
	ca2cb->setBendTowards(hBackbone);
	ca2cb->activate(_sidechain, inherit);

	BondPtr cb2cg1 = BondPtr(new Bond(cBeta, cGamma1));
	cb2cg1->activate(_sidechain, inherit);

	BondPtr cg1hg11 = BondPtr(new Bond(cGamma1, hGamma11));
	cg1hg11->activate(_sidechain, inherit);
	BondPtr cg1hg12 = BondPtr(new Bond(cGamma1, hGamma12));
	cg1hg12->activate(_sidechain, inherit);
	BondPtr cg1hg13 = BondPtr(new Bond(cGamma1, hGamma13));
	cg1hg13->activate(_sidechain, inherit);

	BondPtr cb2cg2 = BondPtr(new Bond(cBeta, cGamma2));
	cb2cg2->activate(_sidechain, inherit);
	BondPtr cb2hb = BondPtr(new Bond(cBeta, hBeta));
	cb2hb->activate(_sidechain, inherit);

	BondPtr cg2hg21 = BondPtr(new Bond(cGamma2, hGamma21));
	cg2hg21->activate(_sidechain, inherit);
	BondPtr cg2hg22 = BondPtr(new Bond(cGamma2, hGamma22));
	cg2hg22->activate(_sidechain, inherit);
	BondPtr cg2hg23 = BondPtr(new Bond(cGamma2, hGamma23));
	cg2hg23->activate(_sidechain, inherit);
}


void Knotter::makeThreonine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr nSpine = backbone->findAtom("N");
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr hBeta = _sidechain->findAtom("HB");
	AtomPtr oGamma1 = _sidechain->findAtom("OG1");
	AtomPtr hGamma11 = _sidechain->findAtom("HG1");
	AtomPtr cGamma2 = _sidechain->findAtom("CG2");
	AtomPtr hGamma21 = _sidechain->findAtom("HG21");
	AtomPtr hGamma22 = _sidechain->findAtom("HG22");
	AtomPtr hGamma23 = _sidechain->findAtom("HG23");
	AtomPtr hBackbone = _sidechain->findAtom("HA");

	AbsolutePtr inherit = std::static_pointer_cast<Absolute>(cAlpha->getModel());

	BondPtr ca2cb = BondPtr(new Bond(cAlpha, cBeta));
	ca2cb->setTorsionAtoms(nSpine, oGamma1);
	ca2cb->setBendTowards(hBackbone);
	ca2cb->activate(_sidechain, inherit);

	BondPtr cb2cg1 = BondPtr(new Bond(cBeta, oGamma1));
	cb2cg1->activate(_sidechain, inherit);
	BondPtr cb2hb = BondPtr(new Bond(cBeta, hBeta));
	cb2hb->activate(_sidechain, inherit);
	BondPtr cb2cg2 = BondPtr(new Bond(cBeta, cGamma2));
	cb2cg2->activate(_sidechain, inherit);

	BondPtr og1hg11 = BondPtr(new Bond(oGamma1, hGamma11));
	og1hg11->activate(_sidechain, inherit);

	BondPtr cg2hg21 = BondPtr(new Bond(cGamma2, hGamma21));
	cg2hg21->activate(_sidechain, inherit);
	BondPtr cg2hg22 = BondPtr(new Bond(cGamma2, hGamma22));
	cg2hg22->activate(_sidechain, inherit);
	BondPtr cg2hg23 = BondPtr(new Bond(cGamma2, hGamma23));
	cg2hg23->activate(_sidechain, inherit);
	
}