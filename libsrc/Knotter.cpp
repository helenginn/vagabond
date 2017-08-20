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
#include "Polymer.h"

#define DEFAULT_PEPTIDE_BLUR 0
#define DEFAULT_RAMACHANDRAN_BLUR 0

void Knotter::tieTowardsCTerminus()
{
	if (!_backbone)
	{
		shout_at_helen("Helen set up a Knotter to tie up atoms,\n"\
					   "and then proceeded to try to convert them\n"\
					   "without specifying a backbone first.\n"\
					   "What an idiot. Haha!");
	}



	PolymerPtr polymer = _backbone->getPolymer();
	int resNum = _backbone->getResNum();
	BackbonePtr prevBackbone, nextBackbone;

	if (resNum > 1)
	{
		MonomerPtr monomer = polymer->getMonomer(resNum - 2);

		if (monomer)
		{
		prevBackbone = monomer->getBackbone(); // one before
		}
	}

	if (resNum < polymer->monomerCount())
	{
		MonomerPtr monomer = polymer->getMonomer(resNum);

		if (monomer)
		{
			nextBackbone = monomer->getBackbone(); // one before
		}
	}

	AtomPtr prevCarbonylCarbon;
	if (prevBackbone)
	{
		prevCarbonylCarbon = prevBackbone->findAtom("C");
	}

	AtomPtr nSpine = _backbone->findAtom("N");
	AtomPtr cAlpha = _backbone->findAtom("CA");
	AtomPtr hAlpha = _backbone->findAtom("HA");

	if (_backbone->getMonomer()->getIdentifier() == "gly")
	{
		hAlpha = _backbone->findAtom("HA3");
	}

	if (_backbone->getMonomer()->getIdentifier() == "pro")
	{
		return;
	}

	AtomPtr carbonylCarbon = _backbone->findAtom("C");
	AtomPtr carbonylOxygen = _backbone->findAtom("O");
	AtomPtr nHydrogen = _backbone->findAtom("H");
	AtomPtr nextNSpine, nextCalpha;

	if (nextBackbone)
	{
		nextNSpine = nextBackbone->findAtom("N");
		nextCalpha = nextBackbone->findAtom("CA");
	}

	AtomPtr inherit = nSpine;

	BondPtr nSpine2cAlpha = BondPtr(new Bond(nSpine, cAlpha));

	if (prevCarbonylCarbon)
	{
		nSpine2cAlpha->setTorsionAtoms(prevCarbonylCarbon, carbonylCarbon);
	}

	nSpine2cAlpha->setBendTowards(nHydrogen);
	nSpine2cAlpha->activate(_backbone, inherit);
	nSpine2cAlpha->addExtraTorsionSample(carbonylOxygen, 0);

	BondPtr cAlpha2Carbonyl = BondPtr(new Bond(cAlpha, carbonylCarbon));
	cAlpha2Carbonyl->setTorsionAtoms(nSpine, carbonylOxygen);
	cAlpha2Carbonyl->activate(_backbone, inherit);

	BondPtr cAlpha2hAlpha = BondPtr(new Bond(cAlpha, hAlpha));
	cAlpha2hAlpha->activate(_backbone, inherit);

	BondPtr carbonyl2oxy = BondPtr(new Bond(carbonylCarbon, carbonylOxygen));
	carbonyl2oxy->activate(_backbone, inherit);

	if (nextBackbone)
	{
		BondPtr carbonyl2nextN = BondPtr(new Bond(carbonylCarbon, nextNSpine));
		carbonyl2nextN->setTorsionAtoms(cAlpha, nextCalpha);
/*
		if (nextCalpha)
		{
			cAlpha2Carbonyl->addExtraTorsionSample(nextCalpha, 0);
		}

		if (nextNSpine)
		{
			nSpine2cAlpha->addExtraTorsionSample(nextNSpine, 0);
		}
*/
		carbonyl2nextN->activate(_backbone, inherit);
	}


}

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

	if (residue == "met")
	{
		convertable = true;
		makeMethionine();
	}

	if (residue == "phe")
	{
		convertable = true;
		makePhenylalanine();
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

	if (residue == "his")
	{
		convertable = true;
		makeHistidine();
	}

	if (residue == "gly")
	{
		convertable = true;
		makeGlycine();
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

void Knotter::makeGlycine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr nSpine = backbone->findAtom("N");
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr hAlpha2 = _sidechain->findAtom("HA2");
	AtomPtr hAlpha3 = backbone->findAtom("HA3");

	AtomPtr inherit = cAlpha;

	BondPtr ca2ha2 = BondPtr(new Bond(cAlpha, hAlpha2));
	ca2ha2->activate(_sidechain, inherit);

}

void Knotter::makeMethionine()
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
	AtomPtr sDelta = _sidechain->findAtom("SD");
	AtomPtr cEpsilon = _sidechain->findAtom("CE");
	AtomPtr hEpsilon1 = _sidechain->findAtom("HE1");
	AtomPtr hEpsilon2 = _sidechain->findAtom("HE2");
	AtomPtr hEpsilon3 = _sidechain->findAtom("HE3");

	AtomPtr inherit = cAlpha;

	BondPtr ca2cb = BondPtr(new Bond(cAlpha, cBeta));
	ca2cb->setBendTowards(hBackbone);
	ca2cb->setTorsionAtoms(nSpine, cGamma);
	ca2cb->activate(_sidechain, inherit);

	BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
	cb2cg->setTorsionAtoms(cAlpha, sDelta);
	cb2cg->activate(_sidechain, inherit);

	BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
	BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));

	BondPtr cg2sd = BondPtr(new Bond(cGamma, sDelta));
	cg2sd->setTorsionAtoms(cBeta, cEpsilon);
	cg2sd->activate(_sidechain, inherit);

	BondPtr cg2hg2 = BondPtr(new Bond(cGamma, hGamma2));
	BondPtr cg2hg3 = BondPtr(new Bond(cGamma, hGamma3));

	BondPtr sd2ce = BondPtr(new Bond(sDelta, cEpsilon));
	sd2ce->activate(_sidechain, inherit);

	BondPtr ce2he1 = BondPtr(new Bond(cEpsilon, hEpsilon1));
	BondPtr ce2he2 = BondPtr(new Bond(cEpsilon, hEpsilon2));
	BondPtr ce2he3 = BondPtr(new Bond(cEpsilon, hEpsilon3));

	ce2he1->activate(_sidechain, inherit);
	ce2he2->activate(_sidechain, inherit);
	ce2he3->activate(_sidechain, inherit);

	cg2hg2->activate(_sidechain, inherit);
	cg2hg3->activate(_sidechain, inherit);

	cb2hb3->activate(_sidechain, inherit);
	cb2hb2->activate(_sidechain, inherit);
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

	AtomPtr inherit = cAlpha;

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

	AtomPtr inherit = cAlpha;

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

	AtomPtr inherit = cAlpha;

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

	AtomPtr inherit = cAlpha;

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

void Knotter::makeHistidine()
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
	AtomPtr cGamma = _sidechain->findAtom("CG");
	AtomPtr nDelta1 = _sidechain->findAtom("ND1");
	AtomPtr hDelta1 = _sidechain->findAtom("HD1");
	AtomPtr cEpsilon1 = _sidechain->findAtom("CE1");
	AtomPtr hEpsilon1 = _sidechain->findAtom("HE1");
	AtomPtr nEpsilon2 = _sidechain->findAtom("NE2");
	AtomPtr hEpsilon2 = _sidechain->findAtom("HE2");
	AtomPtr cDelta2 = _sidechain->findAtom("CD2");
	AtomPtr hDelta2 = _sidechain->findAtom("HD2");
	AtomPtr hBackbone = _sidechain->findAtom("HA");

	AtomPtr inherit = cAlpha;


	BondPtr ca2cb = BondPtr(new Bond(cAlpha, cBeta));
	ca2cb->setBendTowards(hBackbone);
	ca2cb->setTorsionAtoms(nSpine, cGamma);
	ca2cb->activate(_sidechain, inherit);
	ca2cb->addExtraTorsionSample(nDelta1, 0);
	ca2cb->addExtraTorsionSample(cDelta2, 0);
	ca2cb->addExtraTorsionSample(cEpsilon1, 0);
	ca2cb->addExtraTorsionSample(nEpsilon2, 0);

	BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
	cb2cg->setTorsionAtoms(cAlpha, nDelta1);
	cb2cg->activate(_sidechain, inherit);
	cb2cg->addExtraTorsionSample(nEpsilon2, 0);
	cb2cg->addExtraTorsionSample(cEpsilon1, 0);
	cb2cg->addExtraTorsionSample(cDelta2, 0);

	BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
	cb2hb2->activate(_sidechain, inherit);
	BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));
	cb2hb3->activate(_sidechain, inherit);

	BondPtr cg2nd1 = BondPtr(new Bond(cGamma, nDelta1));
	cg2nd1->setFixed(true);
	cg2nd1->setTorsionAtoms(cEpsilon1, cBeta);
	cg2nd1->activate(_sidechain, inherit);

	BondPtr cg2cd2 = BondPtr(new Bond(cGamma, cDelta2));
	cg2cd2->setFixed(true);
	cg2cd2->activate(_sidechain, inherit);

	BondPtr cd22hd2 = BondPtr(new Bond(cDelta2, hDelta2));
	cd22hd2->activate(_sidechain, inherit);

	BondPtr nd12ce1 = BondPtr(new Bond(nDelta1, cEpsilon1));
	nd12ce1->setFixed(true);
	nd12ce1->activate(_sidechain, inherit);
	BondPtr nd12hd1 = BondPtr(new Bond(nDelta1, hDelta1));
	nd12hd1->activate(_sidechain, inherit);


	BondPtr ce12he1 = BondPtr(new Bond(cEpsilon1, hEpsilon1));
	ce12he1->activate(_sidechain, inherit);

	BondPtr cd22ne2 = BondPtr(new Bond(cDelta2, nEpsilon2));
	cd22ne2->setFixed(true);
	cd22ne2->activate(_sidechain, inherit);

	BondPtr ne22he2 = BondPtr(new Bond(nEpsilon2, hEpsilon2));
	ne22he2->activate(_sidechain, inherit);
}


void Knotter::makePhenylalanine()
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
	AtomPtr cGamma = _sidechain->findAtom("CG");
	AtomPtr cDelta1 = _sidechain->findAtom("CD1");
	AtomPtr hDelta1 = _sidechain->findAtom("HD1");
	AtomPtr cEpsilon1 = _sidechain->findAtom("CE1");
	AtomPtr hEpsilon1 = _sidechain->findAtom("HE1");
	AtomPtr cOmega = _sidechain->findAtom("CZ");
	AtomPtr hOmega = _sidechain->findAtom("HZ");
	AtomPtr cEpsilon2 = _sidechain->findAtom("CE2");
	AtomPtr hEpsilon2 = _sidechain->findAtom("HE2");
	AtomPtr cDelta2 = _sidechain->findAtom("CD2");
	AtomPtr hDelta2 = _sidechain->findAtom("HD2");
	AtomPtr hBackbone = _sidechain->findAtom("HA");

	AtomPtr inherit = cAlpha;


	BondPtr ca2cb = BondPtr(new Bond(cAlpha, cBeta));
	ca2cb->setBendTowards(hBackbone);
	ca2cb->setTorsionAtoms(nSpine, cGamma);
	ca2cb->activate(_sidechain, inherit);
	ca2cb->addExtraTorsionSample(cOmega, 0);

	BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
	cb2cg->setTorsionAtoms(cAlpha, cDelta1);
	cb2cg->activate(_sidechain, inherit);
	cb2cg->addExtraTorsionSample(cEpsilon1, 0);
	cb2cg->addExtraTorsionSample(cEpsilon2, 0);

	BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
	cb2hb2->activate(_sidechain, inherit);
	BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));
	cb2hb3->activate(_sidechain, inherit);

	BondPtr cg2cd1 = BondPtr(new Bond(cGamma, cDelta1));
	cg2cd1->setTorsionAtoms(cBeta, cEpsilon1);
	cg2cd1->setFixed(true);
	cg2cd1->activate(_sidechain, inherit);

	BondPtr cg2cd2 = BondPtr(new Bond(cGamma, cDelta2));
	cg2cd2->setTorsionAtoms(cBeta, cEpsilon2);
	cg2cd2->setFixed(true);
	cg2cd2->activate(_sidechain, inherit);

	BondPtr cd22ce2 = BondPtr(new Bond(cDelta2, cEpsilon2));
	cd22ce2->setTorsionAtoms(cGamma, cOmega);
	cd22ce2->setFixed(true);
	cd22ce2->activate(_sidechain, inherit);
	BondPtr cd22hd2 = BondPtr(new Bond(cDelta2, hDelta2));
	cd22hd2->activate(_sidechain, inherit);

	BondPtr cd12ce1 = BondPtr(new Bond(cDelta1, cEpsilon1));
	cd12ce1->setTorsionAtoms(cGamma, hEpsilon1);
	cd12ce1->setFixed(true);
	cd12ce1->activate(_sidechain, inherit);
	BondPtr cd12hd1 = BondPtr(new Bond(cDelta1, hDelta1));
	cd12hd1->activate(_sidechain, inherit);
	BondPtr ce12he1 = BondPtr(new Bond(cEpsilon1, hEpsilon1));
	ce12he1->activate(_sidechain, inherit);

	BondPtr ce22cz = BondPtr(new Bond(cEpsilon2, cOmega));
	ce22cz->setTorsionAtoms(cDelta2, hOmega);
	ce22cz->setFixed(true);
	ce22cz->activate(_sidechain, inherit);

	BondPtr ce22he2 = BondPtr(new Bond(cEpsilon2, hEpsilon2));
	ce22he2->activate(_sidechain, inherit);
	BondPtr cz2hz = BondPtr(new Bond(cOmega, hOmega));

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

	AtomPtr inherit = cAlpha;

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