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
#include "Backbone.h"
#include "Monomer.h"
#include <iostream>
#include "Atom.h"
#include "Bond.h"
#include "GhostBond.h"
#include "Polymer.h"
#include "Options.h"

#define DEFAULT_PEPTIDE_BLUR 0
#define DEFAULT_RAMACHANDRAN_BLUR 0

Knotter::Knotter()
{
	_bondAngles = Options::getBondAngles();
}

BondPtr Knotter::tieBetaCarbon(AtomPtr torsionAtom)
{
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr hAlpha = _sidechain->findAtom("HA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr spineAtom = _backbone->betaCarbonTorsionAtom();

	bool isGlycine = (_backbone->getMonomer()->getIdentifier() == "gly");

	if (isGlycine)
	{
		cBeta= _sidechain->findAtom("HA2");
		hAlpha = _backbone->findAtom("HA3");
	}

	BondPtr ca2cb;

	AtomPtr betaTorsion = _backbone->betaCarbonTorsionAtom();

	if (!betaTorsion)
	{
		shout_at_user("Chain " + _backbone->getPolymer()->getChainID() + ", "
		              "residue " + i_to_str(_backbone->getResNum()) + " is missing the backbone!\n"\
		"Please rebuild and rerun.");
	}

	if (betaTorsion->getAtomName() == "N")
	{
		// tie to C terminus so hAlpha then cBeta
		BondPtr ca2ha = BondPtr(new Bond(cAlpha, hAlpha));
		ca2ha->activate();
		ca2cb = BondPtr(new Bond(cAlpha, cBeta));

		ca2cb->activate();
	}
	else
	{
		// tie to N terminus so cBeta then hAlpha
		ca2cb = BondPtr(new Bond(cAlpha, cBeta));
		ca2cb->activate();
		BondPtr ca2ha = BondPtr(new Bond(cAlpha, hAlpha));
		ca2ha->activate();
	}

	if (_bondAngles > 0)
	{
		ca2cb->setRefineBondAngle();
	}

	return ca2cb;
}

void Knotter::tieTowardsNTerminus()
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

	MonomerPtr monomer = polymer->getMonomer(resNum + 1);

	if (monomer)
	{
		prevBackbone = monomer->getBackbone(); // one before
	}

	if (resNum <= polymer->monomerEnd())
	{
		MonomerPtr monomer = polymer->getMonomer(resNum - 1);

		if (monomer)
		{
			nextBackbone = monomer->getBackbone(); // one before
		}
	}

	AtomPtr prevNitrogen, prevCAlpha;
	if (prevBackbone)
	{
		prevCAlpha = prevBackbone->findAtom("CA");
		prevNitrogen = prevBackbone->findAtom("N");
	}

	AtomPtr nSpine = _backbone->findAtom("N");
	AtomPtr cAlpha = _backbone->findAtom("CA");
	AtomPtr hAlpha = _backbone->findAtom("HA");

	bool isGlycine = (_backbone->getMonomer()->getIdentifier() == "gly");
	if (isGlycine)
	{
		hAlpha = _backbone->findAtom("HA3");
	}

	bool isProline = (_backbone->getMonomer()->getIdentifier() == "pro");
	AtomPtr carbonylCarbon = _backbone->findAtom("C");

	/* already tied up to C terminus */
	if (prevCAlpha->getModel()->isBond())
	{
		BondPtr nitro2Carbon = BondPtr(new Bond(prevNitrogen, carbonylCarbon));
		nitro2Carbon->setRefineFlexibility(false);
		nitro2Carbon->activate();
	}

	AtomPtr carbonylOxygen = _backbone->findAtom("O");
	AtomPtr nHydrogen = _backbone->findAtom("H");

	AtomPtr nextCalpha;

	if (nextBackbone)
	{
		nextCalpha = nextBackbone->findAtom("CA");
	}

	BondPtr carbonyl2CAlpha = BondPtr(new Bond(carbonylCarbon, cAlpha));
	carbonyl2CAlpha->activate();
	
	if (!Options::getPeptideMovement())
	{
		carbonyl2CAlpha->setFixed(true);
		Bond::setTorsion(&*carbonyl2CAlpha, deg2rad(-180));
	}

	BondPtr carbonyl2oxy = BondPtr(new Bond(carbonylCarbon, carbonylOxygen));
	carbonyl2oxy->activate();

	BondPtr cAlpha2NSpine = BondPtr(new Bond(cAlpha, nSpine));

	cAlpha2NSpine->activate();
	
	if (isGlycine && _bondAngles >= 3)
	{
		carbonyl2CAlpha->setRefineBondAngle();
		cAlpha2NSpine->setRefineBondAngle();
	}

	if (nSpine && nSpine->getModel()->isBond())
	{
		BondPtr nSpine2hydrogen = BondPtr(new Bond(nSpine, nHydrogen));
		nSpine2hydrogen->activate();
	}
	
}

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

	if (resNum >= polymer->monomerBegin())
	{
		MonomerPtr monomer = polymer->getMonomer(resNum - 1);

		if (monomer)
		{
			prevBackbone = monomer->getBackbone(); // one before
		}
	}

	if (resNum <= polymer->monomerEnd())
	{
		MonomerPtr monomer = polymer->getMonomer(resNum + 1);

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

	AtomPtr carbonylCarbon = _backbone->findAtom("C");
	AtomPtr carbonylOxygen = _backbone->findAtom("O");
	AtomPtr finalOxygen = _backbone->findAtom("OXT");
	AtomPtr nHydrogen = _backbone->findAtom("H");

	AtomPtr nextNSpine, nextCalpha;

	if (nextBackbone)
	{
		nextNSpine = nextBackbone->findAtom("N");
		nextCalpha = nextBackbone->findAtom("CA");
	}

	BondPtr nSpine2cAlpha = BondPtr(new Bond(nSpine, cAlpha));
	nSpine2cAlpha->activate();
	nSpine2cAlpha->addExtraTorsionSample(carbonylOxygen);

	if (!Options::getPeptideMovement())
	{
		nSpine2cAlpha->setFixed(true);
		Bond::setTorsion(&*nSpine2cAlpha, deg2rad(-180));
	}

	BondPtr cAlpha2Carbonyl = BondPtr(new Bond(cAlpha, carbonylCarbon));
	cAlpha2Carbonyl->activate();

	if (nSpine && nSpine->getModel()->isBond())
	{
		BondPtr nSpine2hydrogen = BondPtr(new Bond(nSpine, nHydrogen));
		nSpine2hydrogen->activate();
	}

	if (nextBackbone)
	{
		BondPtr carbonyl2nextN = BondPtr(new Bond(carbonylCarbon, nextNSpine));
		carbonyl2nextN->setRefineFlexibility(false);
		carbonyl2nextN->activate();
	}

	BondPtr carbonyl2oxy = BondPtr(new Bond(carbonylCarbon, carbonylOxygen));
	carbonyl2oxy->activate();
	
	if (finalOxygen)
	{
		BondPtr carbonyl2oxt = BondPtr(new Bond(carbonylCarbon, finalOxygen));
		carbonyl2oxt->activate();
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

	bool convertible = false;

	if (residue == "ala")
	{
		convertible = true;
		makeAlanine();
	}

	if (residue == "lys")
	{
		convertible = true;
		makeLysine();
	}

	if (residue == "pro")
	{
		convertible = true;
		makeProline();
	}

	if (residue == "ser")
	{
		convertible = true;
		makeSerine();
	}

	if (residue == "met")
	{
		convertible = true;
		makeMethionine();
	}

	if (residue == "leu")
	{
		convertible = true;
		makeLeucine();
	}

	if (residue == "ile")
	{
		convertible = true;
		makeIsoleucine();
	}

	if (residue == "arg")
	{
		convertible = true;
		makeArginine();
	}

	if (residue == "phe")
	{
		convertible = true;
		makePhenylalanine();
	}

	if (residue == "gln")
	{
		convertible = true;
		makeGlutamine();
	}

	if (residue == "glu")
	{
		convertible = true;
		makeGlutamate();
	}

	if (residue == "asn")
	{
		convertible = true;
		makeAsparagine();
	}

	if (residue == "asp")
	{
		convertible = true;
		makeAspartate();
	}
	if (residue == "tyr")
	{
		convertible = true;
		makeTyrosine();
	}

	if (residue == "cys")
	{
		convertible = true;
		makeCysteine();
	}

	if (residue == "val")
	{
		convertible = true;
		makeValine();
	}

	if (residue == "thr")
	{
		convertible = true;
		makeThreonine();
	}

	if (residue == "his")
	{
		convertible = true;
		makeHistidine();
	}

	if (residue == "trp")
	{
		convertible = true;
		makeTryptophan();
	}

	if (residue == "gly")
	{
		convertible = true;
		makeGlycine();
	}

	if (!convertible)
	{
		warn_user("Knotter doesn't know how to tie up " + residue);
	}
}

void Knotter::makeGlycine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	tieBetaCarbon(AtomPtr());
}

void Knotter::makeMethionine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr nAtom = backbone->findAtom("N");
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr hBackbone = _sidechain->findAtom("HA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr cGamma = _sidechain->findAtom("CG");
	AtomPtr hGamma2 = _sidechain->findAtom("HG2");
	AtomPtr hGamma3 = _sidechain->findAtom("HG3");
	AtomPtr sDelta = _sidechain->findAtom("SD");
	AtomPtr cEpsilon = _sidechain->findAtom("CE");
	AtomPtr hEpsilon1 = _sidechain->findAtom("HE1");
	AtomPtr hEpsilon2 = _sidechain->findAtom("HE2");
	AtomPtr hEpsilon3 = _sidechain->findAtom("HE3");

	tieBetaCarbon(cGamma);

	BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
	cb2cg->activate();

	BondPtr cg2sd = BondPtr(new Bond(cGamma, sDelta));
	cg2sd->activate();

	if (_bondAngles >= 2)
	{
		cb2cg->setRefineBondAngle();
		cg2sd->setRefineBondAngle();
	}

	BondPtr sd2ce = BondPtr(new Bond(sDelta, cEpsilon));
	sd2ce->activate();

}

void Knotter::makeArginine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr nAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr hBackbone = _sidechain->findAtom("HA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr cGamma = _sidechain->findAtom("CG");
	AtomPtr cDelta = _sidechain->findAtom("CD");
	AtomPtr nEpsilon = _sidechain->findAtom("NE");
	AtomPtr hEpsilon = _sidechain->findAtom("HE");
	AtomPtr cOmega = _sidechain->findAtom("CZ");
	AtomPtr nOmega1 = _sidechain->findAtom("NH1");
	AtomPtr nOmega2 = _sidechain->findAtom("NH2");

	tieBetaCarbon(cGamma);

	BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
	cb2cg->activate();

	BondPtr cg2cd = BondPtr(new Bond(cGamma, cDelta));
	cg2cd->activate();

	BondPtr cd2ce = BondPtr(new Bond(cDelta, nEpsilon));
	cd2ce->activate();

	BondPtr ne2cz = BondPtr(new Bond(nEpsilon, cOmega));
	ne2cz->activate();

	BondPtr ce2nh1= BondPtr(new Bond(cOmega, nOmega1));
	BondPtr ce2nh2 = BondPtr(new Bond(cOmega, nOmega2));
	ce2nh1->activate();
	ce2nh2->activate();
}

void Knotter::makeLysine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr hBackbone = _sidechain->findAtom("HA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr cGamma = _sidechain->findAtom("CG");
	AtomPtr cDelta = _sidechain->findAtom("CD");
	AtomPtr cEpsilon = _sidechain->findAtom("CE");
	AtomPtr nOmega = _sidechain->findAtom("NZ");

	tieBetaCarbon(cGamma);

	BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
	cb2cg->activate();

	BondPtr cg2cd = BondPtr(new Bond(cGamma, cDelta));
	cg2cd->activate();

	BondPtr cd2ce = BondPtr(new Bond(cDelta, cEpsilon));
	cd2ce->activate();

	BondPtr ce2nz = BondPtr(new Bond(cEpsilon, nOmega));
	ce2nz->activate();
}

void Knotter::makeProline()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr nAtom = backbone->findAtom("N");
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr hBackbone = _sidechain->findAtom("HA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr cGamma = _sidechain->findAtom("CG");
	AtomPtr cDelta = _sidechain->findAtom("CD");

	BondPtr ca2cb = tieBetaCarbon(cGamma);
	ca2cb->setRefineBondAngle(false);
	ca2cb->setRefineFlexibility(false);

	BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
	cb2cg->setRefineFlexibility(false);
	ca2cb->addExtraTorsionSample(cGamma);
	cb2cg->activate();

	BondPtr cg2cd = BondPtr(new Bond(cGamma, cDelta));
	cg2cd->setRefineFlexibility(false);
	ca2cb->addExtraTorsionSample(cDelta);
	cg2cd->activate();
	
	GhostBondPtr ghost = GhostBondPtr(new GhostBond());
	ghost->setAtoms(cDelta, nAtom);
}


void Knotter::makeSerine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);
	//int resNum = monomer->getResidueNum();

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr oGamma = _sidechain->findAtom("OG");
	AtomPtr hGamma = _sidechain->findAtom("HG");

	tieBetaCarbon(oGamma);

	BondPtr cb2og = BondPtr(new Bond(cBeta, oGamma));
	cb2og->activate();

	BondPtr og2hg = BondPtr(new Bond(oGamma, hGamma));
	og2hg->activate();
}

void Knotter::makeCysteine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr sGamma = _sidechain->findAtom("SG");
	AtomPtr hGamma = _sidechain->findAtom("HG");

	tieBetaCarbon(sGamma);

	BondPtr cb2sg = BondPtr(new Bond(cBeta, sGamma));
	cb2sg->activate();
	if (_bondAngles >= 2)
	{
		cb2sg->setRefineBondAngle();
	}

	BondPtr sg2hg = BondPtr(new Bond(sGamma, hGamma));
	sg2hg->activate();
}

void Knotter::makeValine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr cGamma1 = _sidechain->findAtom("CG1");
	AtomPtr cGamma2 = _sidechain->findAtom("CG2");

	tieBetaCarbon(cGamma1);

	BondPtr cb2cg1 = BondPtr(new Bond(cBeta, cGamma1));
	cb2cg1->activate();

	BondPtr cb2cg2 = BondPtr(new Bond(cBeta, cGamma2));
	cb2cg2->activate();
}

void Knotter::makeAlanine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");

	tieBetaCarbon(AtomPtr());
}


void Knotter::makeHistidine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr cGamma = _sidechain->findAtom("CG");
	AtomPtr nDelta1 = _sidechain->findAtom("ND1");
	AtomPtr cEpsilon1 = _sidechain->findAtom("CE1");
	AtomPtr nEpsilon2 = _sidechain->findAtom("NE2");
	AtomPtr cDelta2 = _sidechain->findAtom("CD2");

	BondPtr ca2cb = tieBetaCarbon(cGamma);

	BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
	cb2cg->activate();

	if (_bondAngles >= 1)
	{
		cb2cg->setRefineBondAngle();
	}

	BondPtr cg2nd1 = BondPtr(new Bond(cGamma, nDelta1));
	cg2nd1->activate();
	cg2nd1->addExtraTorsionSample(cDelta2);
	cg2nd1->addExtraTorsionSample(cEpsilon1);
	cg2nd1->addExtraTorsionSample(nEpsilon2);

	BondPtr cg2cd2 = BondPtr(new Bond(cGamma, cDelta2));
	cg2cd2->setFixed(true);
	cg2cd2->activate();

	BondPtr nd12ce1 = BondPtr(new Bond(nDelta1, cEpsilon1));
	nd12ce1->setFixed(true);
	nd12ce1->activate();

	BondPtr cd22ne2 = BondPtr(new Bond(cDelta2, nEpsilon2));
	cd22ne2->setFixed(true);
	cd22ne2->activate();

	GhostBondPtr ghost = GhostBondPtr(new GhostBond());
	ghost->setAtoms(nEpsilon2, cEpsilon1);
}

void Knotter::makeTyrosine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr cGamma = _sidechain->findAtom("CG");
	AtomPtr cDelta1 = _sidechain->findAtom("CD1");
	AtomPtr cEpsilon1 = _sidechain->findAtom("CE1");
	AtomPtr cOmega = _sidechain->findAtom("CZ");
	AtomPtr oxygen = _sidechain->findAtom("OH");
	AtomPtr cEpsilon2 = _sidechain->findAtom("CE2");
	AtomPtr cDelta2 = _sidechain->findAtom("CD2");

	BondPtr ca2cb = tieBetaCarbon(cGamma);

	BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
	cb2cg->activate();

	if (_bondAngles >= 1)
	{
		ca2cb->setRefineBondAngle();
		cb2cg->setRefineBondAngle();
	}

	BondPtr cg2cd1 = BondPtr(new Bond(cGamma, cDelta1));
	cg2cd1->activate();
	cg2cd1->addExtraTorsionSample(cEpsilon1);
	cg2cd1->addExtraTorsionSample(cEpsilon2);
	cg2cd1->addExtraTorsionSample(cDelta2);
	cg2cd1->addExtraTorsionSample(cOmega);
	cg2cd1->addExtraTorsionSample(oxygen);

	BondPtr cg2cd2 = BondPtr(new Bond(cGamma, cDelta2));
	cg2cd2->setFixed(true);
	cg2cd2->activate();

	BondPtr cd22ce2 = BondPtr(new Bond(cDelta2, cEpsilon2));
	cd22ce2->setFixed(true);
	cd22ce2->activate();

	BondPtr cd12ce1 = BondPtr(new Bond(cDelta1, cEpsilon1));
	cd12ce1->setFixed(true);
	cd12ce1->activate();
	BondPtr ce22cz = BondPtr(new Bond(cEpsilon2, cOmega));
	ce22cz->setFixed(true);
	ce22cz->activate();

	BondPtr cz2oh = BondPtr(new Bond(cOmega, oxygen));
	cz2oh->activate();

	GhostBondPtr ghost = GhostBondPtr(new GhostBond());
	ghost->setAtoms(cEpsilon1, cOmega);
}

void Knotter::makePhenylalanine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr cGamma = _sidechain->findAtom("CG");
	AtomPtr cDelta1 = _sidechain->findAtom("CD1");
	AtomPtr cEpsilon1 = _sidechain->findAtom("CE1");
	AtomPtr cOmega = _sidechain->findAtom("CZ");
	AtomPtr hOmega = _sidechain->findAtom("HZ");
	AtomPtr cEpsilon2 = _sidechain->findAtom("CE2");
	AtomPtr cDelta2 = _sidechain->findAtom("CD2");

	BondPtr ca2cb = tieBetaCarbon(cGamma);

	BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
	cb2cg->addExtraTorsionSample(cOmega);
	cb2cg->activate();


	if (_bondAngles >= 1)
	{
		ca2cb->setRefineBondAngle();
		cb2cg->setRefineBondAngle();
	}

	BondPtr cg2cd1 = BondPtr(new Bond(cGamma, cDelta1));
	cg2cd1->activate();
	cg2cd1->addExtraTorsionSample(cEpsilon1);
	cg2cd1->addExtraTorsionSample(cEpsilon2);
	cg2cd1->addExtraTorsionSample(cDelta2);
	cg2cd1->addExtraTorsionSample(cOmega);

	BondPtr cg2cd2 = BondPtr(new Bond(cGamma, cDelta2));
	cg2cd2->setFixed(true);
	cg2cd2->activate();

	BondPtr cd22ce2 = BondPtr(new Bond(cDelta2, cEpsilon2));
	cd22ce2->setFixed(true);
	cd22ce2->activate();

	BondPtr cd12ce1 = BondPtr(new Bond(cDelta1, cEpsilon1));
	cd12ce1->setFixed(true);
	cd12ce1->activate();

	BondPtr ce22cz = BondPtr(new Bond(cEpsilon2, cOmega));
	ce22cz->setFixed(true);
	ce22cz->activate();

	BondPtr cz2hz = BondPtr(new Bond(cOmega, hOmega));
	cz2hz->activate();

	GhostBondPtr ghost = GhostBondPtr(new GhostBond());
	ghost->setAtoms(cEpsilon1, cOmega);
}

void Knotter::makeTryptophan()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr cGamma = _sidechain->findAtom("CG");
	AtomPtr cDelta1 = _sidechain->findAtom("CD1");
	AtomPtr cDelta2 = _sidechain->findAtom("CD2");
	AtomPtr nEpsilon1 = _sidechain->findAtom("NE1");
	AtomPtr cEpsilon2 = _sidechain->findAtom("CE2");
	AtomPtr cEpsilon3 = _sidechain->findAtom("CE3");
	AtomPtr cOmega2 = _sidechain->findAtom("CZ2");
	AtomPtr cOmega3 = _sidechain->findAtom("CZ3");
	AtomPtr ch2 = _sidechain->findAtom("CH2");

	BondPtr ca2cb = tieBetaCarbon(cGamma);

	BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
	cb2cg->activate();

	if (_bondAngles >= 1)
	{
		ca2cb->setRefineBondAngle();
		cb2cg->setRefineBondAngle();
	}

	BondPtr cg2cd1 = BondPtr(new Bond(cGamma, cDelta1));
	cg2cd1->activate();
	cg2cd1->addExtraTorsionSample(cDelta2);
	cg2cd1->addExtraTorsionSample(nEpsilon1);
	cg2cd1->addExtraTorsionSample(cEpsilon2);
	cg2cd1->addExtraTorsionSample(cEpsilon3);
	cg2cd1->addExtraTorsionSample(cOmega2);
	cg2cd1->addExtraTorsionSample(cOmega3);
	cg2cd1->addExtraTorsionSample(ch2);

	BondPtr cg2cd2 = BondPtr(new Bond(cGamma, cDelta2));
	cg2cd2->setFixed(true);
	cg2cd2->activate();

	BondPtr cd22ce2 = BondPtr(new Bond(cDelta2, cEpsilon3));
	cd22ce2->setFixed(true);
	cd22ce2->activate();

	BondPtr cd12ce1 = BondPtr(new Bond(cDelta1, nEpsilon1));
	cd12ce1->setFixed(true);
	cd12ce1->activate();

	BondPtr ne12ce2 = BondPtr(new Bond(nEpsilon1, cEpsilon2));
	ne12ce2->activate();

	BondPtr ce22cz = BondPtr(new Bond(cEpsilon2, cOmega2));
	ce22cz->setFixed(true);
	ce22cz->activate();

	BondPtr ce22ce2 = BondPtr(new Bond(cEpsilon3, cOmega3));
	ce22ce2->activate();

	BondPtr ce32ce3 = BondPtr(new Bond(cOmega3, ch2));
	ce32ce3->activate();

	GhostBondPtr ghost = GhostBondPtr(new GhostBond());
	ghost->setAtoms(cOmega2, ch2);

	ghost = GhostBondPtr(new GhostBond());
	ghost->setAtoms(cDelta2, cEpsilon2);
}

void Knotter::makeIsoleucine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr cGamma1 = _sidechain->findAtom("CG1");
	AtomPtr cGamma2 = _sidechain->findAtom("CG2");
	AtomPtr cDelta1 = _sidechain->findAtom("CD1");

	BondPtr ca2cb = tieBetaCarbon(cGamma1);

	BondPtr cb2cg1 = BondPtr(new Bond(cBeta, cGamma1));
	cb2cg1->activate();
	BondPtr cb2cg2 = BondPtr(new Bond(cBeta, cGamma2));
	cb2cg2->activate();

	BondPtr cg1cd11 = BondPtr(new Bond(cGamma1, cDelta1));
	cg1cd11->activate();
}

void Knotter::makeLeucine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr cGamma = _sidechain->findAtom("CG");
	AtomPtr cDelta1 = _sidechain->findAtom("CD1");
	AtomPtr cDelta2 = _sidechain->findAtom("CD2");

	BondPtr ca2cb = tieBetaCarbon(cGamma);

	BondPtr cb2cg1 = BondPtr(new Bond(cBeta, cGamma));
	cb2cg1->activate();
	BondPtr cg1cd1 = BondPtr(new Bond(cGamma, cDelta1));
	cg1cd1->activate();
	BondPtr cg1cd2 = BondPtr(new Bond(cGamma, cDelta2));
	cg1cd2->activate();
}


void Knotter::makeAspartate()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr cGamma = _sidechain->findAtom("CG");
	AtomPtr oDelta1 = _sidechain->findAtom("OD1");
	AtomPtr oDelta2 = _sidechain->findAtom("OD2");

	BondPtr ca2cb = tieBetaCarbon(cGamma);

	BondPtr cb2cg1 = BondPtr(new Bond(cBeta, cGamma));
	cb2cg1->activate();

	BondPtr cg1cd1 = BondPtr(new Bond(cGamma, oDelta1));
	cg1cd1->activate();
	BondPtr cg1cd2 = BondPtr(new Bond(cGamma, oDelta2));
	cg1cd2->activate();
}

void Knotter::makeAsparagine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr cGamma = _sidechain->findAtom("CG");
	AtomPtr oDelta1 = _sidechain->findAtom("OD1");
	AtomPtr nDelta2 = _sidechain->findAtom("ND2");

	BondPtr ca2cb = tieBetaCarbon(cGamma);

	BondPtr cb2cg1 = BondPtr(new Bond(cBeta, cGamma));
	cb2cg1->activate();

	BondPtr cg1nd2 = BondPtr(new Bond(cGamma, nDelta2));
	cg1nd2->activate();
	BondPtr cg1od1 = BondPtr(new Bond(cGamma, oDelta1));
	cg1od1->activate();
}

void Knotter::makeGlutamine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr cGamma = _sidechain->findAtom("CG");
	AtomPtr cDelta = _sidechain->findAtom("CD");
	AtomPtr oEpsilon1 = _sidechain->findAtom("OE1");
	AtomPtr nEpsilon2 = _sidechain->findAtom("NE2");

	BondPtr ca2cb = tieBetaCarbon(cGamma);

	BondPtr cb2cg1 = BondPtr(new Bond(cBeta, cGamma));
	cb2cg1->activate();

	BondPtr cg2cd1 = BondPtr(new Bond(cGamma, cDelta));
	cg2cd1->activate();

	BondPtr cd2ce1 = BondPtr(new Bond(cDelta, oEpsilon1));
	cd2ce1->activate();
	BondPtr cd2ce2 = BondPtr(new Bond(cDelta, nEpsilon2));
	cd2ce2->activate();
}


void Knotter::makeGlutamate()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr cGamma = _sidechain->findAtom("CG");
	AtomPtr cDelta = _sidechain->findAtom("CD");
	AtomPtr oEpsilon1 = _sidechain->findAtom("OE1");
	AtomPtr oEpsilon2 = _sidechain->findAtom("OE2");

	BondPtr ca2cb = tieBetaCarbon(cGamma);

	BondPtr cb2cg1 = BondPtr(new Bond(cBeta, cGamma));
	cb2cg1->activate();

	BondPtr cg2cd1 = BondPtr(new Bond(cGamma, cDelta));
	cg2cd1->activate();

	BondPtr cd2ce1 = BondPtr(new Bond(cDelta, oEpsilon1));
	cd2ce1->activate();
	BondPtr cd2ce2 = BondPtr(new Bond(cDelta, oEpsilon2));
	cd2ce2->activate();
}

void Knotter::makeThreonine()
{
	MonomerPtr monomer = _sidechain->getMonomer();
	BackbonePtr backbone = monomer->getBackbone();
	std::string residue = monomer->getIdentifier();
	_sidechain->setCanRefine(true);

	AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
	AtomPtr cAlpha = _sidechain->findAtom("CA");
	AtomPtr cBeta = _sidechain->findAtom("CB");
	AtomPtr oGamma1 = _sidechain->findAtom("OG1");
	AtomPtr cGamma2 = _sidechain->findAtom("CG2");

	BondPtr ca2cb = tieBetaCarbon(oGamma1);

	BondPtr cb2cg1 = BondPtr(new Bond(cBeta, oGamma1));
	cb2cg1->activate();
	BondPtr cb2cg2 = BondPtr(new Bond(cBeta, cGamma2));
	cb2cg2->activate();
}
