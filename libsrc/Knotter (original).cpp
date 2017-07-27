//
//  Knotter.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Knotter.h"
#include "Shouter.h"
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

	if (resNum > 120 || resNum == 58)
	{
		return;
	}
	
	if (residue == "ser")
	{
		convertable = true;
		makeSerine();
	}

	if (convertable)
	{
		std::cout << "Knotter is converting residue: " << residue
		<< resNum << std::endl;
	}
	else
	{
		std::cout << "Knotter doesn't know how to tie up " << residue << std::endl;
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
	ca2cb->setTorsionAtoms(nSpine, oGamma);
	ca2cb->activate();

	BondPtr cb2og = BondPtr(new Bond(cBeta, oGamma));
	Bond::setBondLength(&*cb2og, 1.41884);
	cb2og->activate();

	BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
	cb2hb2->activate();
	BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));
	cb2hb3->activate();

	BondPtr og2hg = BondPtr(new Bond(oGamma, hGamma));
	og2hg->activate();






}