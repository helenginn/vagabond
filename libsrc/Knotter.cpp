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

	if (resNum < 23 || resNum > 104)
	{
		return;
	}

	if (residue == "ser")
	{
		convertable = true;
		makeSerine();
	}

	if (residue == "val")
	{
		convertable = true;
		makeValine();
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
	Bond::setBondLength(&*ca2cb, 1.52293);
	Bond::setTorsionBlur(&*ca2cb, deg2rad(1.0));
	ca2cb->setTorsionAtoms(nSpine, oGamma);
	ca2cb->activate(_sidechain, inherit);

	BondPtr cb2og = BondPtr(new Bond(cBeta, oGamma));
	Bond::setBondLength(&*cb2og, 1.41884);
	cb2og->activate(_sidechain, inherit);

	BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
	cb2hb2->activate(_sidechain, inherit);
	BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));
	cb2hb3->activate(_sidechain, inherit);

	BondPtr og2hg = BondPtr(new Bond(oGamma, hGamma));
	og2hg->activate(_sidechain, inherit);
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

	AbsolutePtr inherit = std::static_pointer_cast<Absolute>(cAlpha->getModel());

	if (!hBeta || !hGamma11 || !hGamma12 || !hGamma13 || !hGamma21
		|| !hGamma22 || !hGamma23)
	{
		warn_user("Missing some hydrogens in " + residue + ", continuing.");
	}

	if (!cAlpha || !cBeta || !cGamma1 || !cGamma2)
	{
		warn_user("Missing vital atoms in " + residue + ", giving up.");
		return;
	}

	BondPtr ca2cb = BondPtr(new Bond(cAlpha, cBeta));
	Bond::setBondLength(&*ca2cb, 1.52349);
	ca2cb->setTorsionAtoms(nSpine, cGamma1);
	Bond::setTorsionBlur(&*ca2cb, deg2rad(2.0));
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