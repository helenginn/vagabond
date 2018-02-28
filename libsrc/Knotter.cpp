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
#include "Polymer.h"
#include "Options.h"

#define DEFAULT_PEPTIDE_BLUR 0
#define DEFAULT_RAMACHANDRAN_BLUR 0

BondPtr Knotter::tieBetaCarbon(AtomPtr torsionAtom)
{
    AtomPtr cAlpha = _sidechain->findAtom("CA");
    AtomPtr hAlpha = _sidechain->findAtom("HA");
    AtomPtr cBeta = _sidechain->findAtom("CB");
    AtomPtr spineAtom = _backbone->betaCarbonTorsionAtom();

    if (_backbone->getMonomer()->getIdentifier() == "gly")
    {
        cBeta= _sidechain->findAtom("HA2");
        hAlpha = _backbone->findAtom("HA3");
    }

    BondPtr ca2cb;

    AtomPtr betaTorsion = _backbone->betaCarbonTorsionAtom();

    if (!betaTorsion)
    {
        shout_at_user("Residue " + _backbone->getPolymer()->getChainID()
                      + i_to_str(_backbone->getResNum()) + " is missing the backbone!\n"\
                      "Please rebuild and rerun.");
    }

    if (betaTorsion->getAtomName() == "N")
    {
        // tie to C terminus so hAlpha then cBeta
        BondPtr ca2ha = BondPtr(new Bond(cAlpha, hAlpha));
        ca2ha->activate(_sidechain);
        ca2cb = BondPtr(new Bond(cAlpha, cBeta));
        if (torsionAtom && spineAtom)
        {
            ca2cb->setTorsionAtoms(spineAtom, torsionAtom);
        }
        ca2cb->activate(_sidechain);
    }
    else
    {
        // tie to N terminus so cBeta then hAlpha
        ca2cb = BondPtr(new Bond(cAlpha, cBeta));
        if (torsionAtom && spineAtom)
        {
            ca2cb->setTorsionAtoms(spineAtom, torsionAtom);
        }
        ca2cb->activate(_sidechain);
        BondPtr ca2ha = BondPtr(new Bond(cAlpha, hAlpha));
        ca2ha->activate(_sidechain);
    }

    ca2cb->setRefineBondAngle();

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

//    if (resNum > 1)
    {
        /* to get the previous (C-terminal) backbone residue */
        MonomerPtr monomer = polymer->getMonomer(resNum + 1);

        if (monomer)
        {
            prevBackbone = monomer->getBackbone(); // one before
        }
    }

    if (resNum < polymer->monomerCount())
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

    if (_backbone->getMonomer()->getIdentifier() == "gly")
    {
        hAlpha = _backbone->findAtom("HA3");
    }

    bool isProline = (_backbone->getMonomer()->getIdentifier() == "pro");
    AtomPtr carbonylCarbon = _backbone->findAtom("C");

    if (prevCAlpha->getModel()->getClassName() == "Bond")
    {
        BondPtr nitro2Carbon = BondPtr(new Bond(prevNitrogen, carbonylCarbon));
//		nitro2Carbon->setTorsionStepMult(0.2);
        nitro2Carbon->setTorsionAtoms(prevCAlpha, cAlpha);
        nitro2Carbon->activate(_backbone);
    //    Bond::setTorsion(&*nitro2Carbon, deg2rad(180));
    }

    AtomPtr carbonylOxygen = _backbone->findAtom("O");
    AtomPtr nHydrogen = _backbone->findAtom("H");

    AtomPtr nextCSpine, nextCalpha;

    if (nextBackbone)
    {
        nextCSpine = nextBackbone->findAtom("C");
        nextCalpha = nextBackbone->findAtom("CA");
    }

    AtomPtr inherit = nSpine;

    BondPtr carbonyl2CAlpha = BondPtr(new Bond(carbonylCarbon, cAlpha));
    carbonyl2CAlpha->setTorsionAtoms(prevNitrogen, nSpine);
    carbonyl2CAlpha->activate(_backbone);

    BondPtr carbonyl2oxy = BondPtr(new Bond(carbonylCarbon, carbonylOxygen));
    carbonyl2oxy->activate(_backbone);

    BondPtr cAlpha2NSpine = BondPtr(new Bond(cAlpha, nSpine));

    if (nextCalpha)
    {
        cAlpha2NSpine->setTorsionAtoms(carbonylCarbon, nextCalpha);
    }

    cAlpha2NSpine->activate(_backbone);

    if (nSpine && nSpine->getModel()->isBond())
    {
        BondPtr nSpine2hydrogen = BondPtr(new Bond(nSpine, nHydrogen));
        nSpine2hydrogen->activate(_backbone);
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

    if (resNum > 1)
    {
        MonomerPtr monomer = polymer->getMonomer(resNum - 1);

        if (monomer)
        {
            prevBackbone = monomer->getBackbone(); // one before
        }
    }

    if (resNum < polymer->monomerCount())
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

    nSpine2cAlpha->activate(_backbone, inherit);

    if (nSpine && nSpine->getModel()->isBond())
    {
        BondPtr nSpine2hydrogen = BondPtr(new Bond(nSpine, nHydrogen));
        nSpine2hydrogen->activate(_backbone);
    }


    nSpine2cAlpha->addExtraTorsionSample(carbonylOxygen, 0);

    BondPtr cAlpha2Carbonyl = BondPtr(new Bond(cAlpha, carbonylCarbon));
    if (nextBackbone && nSpine)
    {
        cAlpha2Carbonyl->setTorsionAtoms(nSpine, nextNSpine);
    }
    else
    {
        cAlpha2Carbonyl->setTorsionAtoms(nSpine, carbonylOxygen);
    }

    cAlpha2Carbonyl->activate(_backbone, inherit);

//    BondPtr cAlpha2hAlpha = BondPtr(new Bond(cAlpha, hAlpha));
//    cAlpha2hAlpha->activate(_backbone, inherit);

    if (nextBackbone)
    {
        BondPtr carbonyl2nextN = BondPtr(new Bond(carbonylCarbon, nextNSpine));
        carbonyl2nextN->setTorsionAtoms(cAlpha, nextCalpha);
        carbonyl2nextN->activate(_backbone, inherit);
//		carbonyl2nextN->setTorsionStepMult(0.2);
//        Bond::setTorsion(&*carbonyl2nextN, deg2rad(180));

    }

    BondPtr carbonyl2oxy = BondPtr(new Bond(carbonylCarbon, carbonylOxygen));
    carbonyl2oxy->activate(_backbone, inherit);
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

    bool convertable = false;

    if (residue == "ala")
    {
        convertable = true;
        makeAlanine();
    }

    if (residue == "lys")
    {
        convertable = true;
        makeLysine();
    }

    if (residue == "pro")
    {
        convertable = true;
        makeProline();
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

    if (residue == "leu")
    {
        convertable = true;
        makeLeucine();
    }

    if (residue == "ile")
    {
        convertable = true;
        makeIsoleucine();
    }

    if (residue == "arg")
    {
        convertable = true;
        makeArginine();
    }

    if (residue == "phe")
    {
        convertable = true;
        makePhenylalanine();
    }

    if (residue == "gln")
    {
        convertable = true;
        makeGlutamine();
    }

    if (residue == "glu")
    {
        convertable = true;
        makeGlutamate();
    }

    if (residue == "asn")
    {
        convertable = true;
        makeAsparagine();
    }

    if (residue == "asp")
    {
        convertable = true;
        makeAspartate();
    }
    if (residue == "tyr")
    {
        convertable = true;
        makeTyrosine();
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

    if (residue == "trp")
    {
        convertable = true;
        makeTryptophan();
    }

    if (residue == "gly")
    {
        convertable = true;
        makeGlycine();
    }

    if (!convertable)
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

    tieBetaCarbon(cGamma);

    BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
    cb2cg->setTorsionAtoms(cAlpha, sDelta);
    cb2cg->activate(_sidechain, inherit);

    BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
    BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));

    BondPtr cg2sd = BondPtr(new Bond(cGamma, sDelta));
    cg2sd->setTorsionAtoms(cBeta, cEpsilon);
    cg2sd->activate(_sidechain, inherit);

    cb2cg->setRefineBondAngle();
    cg2sd->setRefineBondAngle();

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

void Knotter::makeArginine()
{
    MonomerPtr monomer = _sidechain->getMonomer();
    BackbonePtr backbone = monomer->getBackbone();
    std::string residue = monomer->getIdentifier();
    _sidechain->setCanRefine(true);

    AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
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
    AtomPtr nEpsilon = _sidechain->findAtom("NE");
    AtomPtr hEpsilon = _sidechain->findAtom("HE");
    AtomPtr cOmega = _sidechain->findAtom("CZ");
    AtomPtr nOmega1 = _sidechain->findAtom("NH1");
    AtomPtr nOmega2 = _sidechain->findAtom("NH2");
    AtomPtr hh11 = _sidechain->findAtom("HH11");
    AtomPtr hh12 = _sidechain->findAtom("HH12");
    AtomPtr hh21 = _sidechain->findAtom("HH21");
    AtomPtr hh22 = _sidechain->findAtom("HH22");

    AtomPtr inherit = cAlpha;

    tieBetaCarbon(cGamma);

    BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
    cb2cg->setTorsionAtoms(cAlpha, cDelta);
    cb2cg->activate(_sidechain, inherit);

    BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
    BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));

    BondPtr cg2cd = BondPtr(new Bond(cGamma, cDelta));
    cg2cd->setTorsionAtoms(cBeta, nEpsilon);
    cg2cd->activate(_sidechain, inherit);

    BondPtr cg2hg2 = BondPtr(new Bond(cGamma, hGamma2));
    BondPtr cg2hg3 = BondPtr(new Bond(cGamma, hGamma3));

    BondPtr cd2ce = BondPtr(new Bond(cDelta, nEpsilon));
    cd2ce->setTorsionAtoms(cGamma, cOmega);
    cd2ce->activate(_sidechain, inherit);

    BondPtr cd2hd2 = BondPtr(new Bond(cDelta, hDelta2));
    BondPtr cd2hd3 = BondPtr(new Bond(cDelta, hDelta3));

    BondPtr ne2cz = BondPtr(new Bond(nEpsilon, cOmega));
    BondPtr ne2he = BondPtr(new Bond(nEpsilon, hEpsilon));
    ne2cz->setTorsionAtoms(cDelta, nOmega1);
    ne2cz->activate(_sidechain, inherit);
    ne2he->activate(_sidechain, inherit);

    BondPtr ce2nh1= BondPtr(new Bond(cOmega, nOmega1));
    BondPtr ce2nh2 = BondPtr(new Bond(cOmega, nOmega2));
    ce2nh1->activate(_sidechain, inherit);
    ce2nh2->activate(_sidechain, inherit);

    BondPtr nz1hz11 = BondPtr(new Bond(nOmega1, hh11));
    BondPtr nz1hz12 = BondPtr(new Bond(nOmega1, hh12));
    BondPtr nz2hz21 = BondPtr(new Bond(nOmega2, hh21));
    BondPtr nz2hz22 = BondPtr(new Bond(nOmega2, hh22));

    nz1hz11->activate(_sidechain, inherit);
    nz1hz11->activate(_sidechain, inherit);
    nz2hz21->activate(_sidechain, inherit);
    nz2hz22->activate(_sidechain, inherit);


    cd2hd2->activate(_sidechain, inherit);
    cd2hd3->activate(_sidechain, inherit);

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

    AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
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

    tieBetaCarbon(cGamma);

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

void Knotter::makeProline()
{
    MonomerPtr monomer = _sidechain->getMonomer();
    BackbonePtr backbone = monomer->getBackbone();
    std::string residue = monomer->getIdentifier();
    _sidechain->setCanRefine(true);

    AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
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

    AtomPtr inherit = cAlpha;

    BondPtr ca2cb = tieBetaCarbon(cGamma);
    ca2cb->addExtraTorsionSample(cGamma, 0);
    ca2cb->addExtraTorsionSample(cDelta, 0);
    double torsion = deg2rad(-92.16);
    if (_backbone->betaCarbonTorsionAtom()->getAtomName() == "N")
    {
        torsion = deg2rad(27.27);
    }
//    Bond::setTorsion(&*ca2cb, torsion);
    ca2cb->setRefineBondAngle(false);

    BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
//    cb2cg->setFixed(true);
    cb2cg->setTorsionAtoms(cAlpha, cDelta);
    cb2cg->activate(_sidechain, inherit);
//    Bond::setTorsion(&*cb2cg, deg2rad(-34.61));

    BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
    BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));
    cb2hb2->activate(_sidechain, inherit);
    cb2hb3->activate(_sidechain, inherit);

    BondPtr cg2cd = BondPtr(new Bond(cGamma, cDelta));
//    cg2cd->setFixed(true);
    cg2cd->activate(_sidechain, inherit);

    BondPtr cg2hg2 = BondPtr(new Bond(cGamma, hGamma2));
    cg2hg2->activate(_sidechain, inherit);
    BondPtr cg2hg3 = BondPtr(new Bond(cGamma, hGamma3));
    cg2hg3->activate(_sidechain, inherit);
    BondPtr cd2hd2 = BondPtr(new Bond(cDelta, hDelta2));
    cd2hd2->activate(_sidechain, inherit);
    BondPtr cd2hd3 = BondPtr(new Bond(cDelta, hDelta3));
    cd2hd3->activate(_sidechain, inherit);
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
    AtomPtr hBeta2 = _sidechain->findAtom("HB2");
    AtomPtr hBeta3 = _sidechain->findAtom("HB3");
    AtomPtr oGamma = _sidechain->findAtom("OG");
    AtomPtr hGamma = _sidechain->findAtom("HG");
    AtomPtr hBackBone = _sidechain->findAtom("HA");

    AtomPtr inherit = cAlpha;

    tieBetaCarbon(oGamma);

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

    AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
    AtomPtr cAlpha = _sidechain->findAtom("CA");
    AtomPtr cBeta = _sidechain->findAtom("CB");
    AtomPtr hBeta2 = _sidechain->findAtom("HB2");
    AtomPtr hBeta3 = _sidechain->findAtom("HB3");
    AtomPtr sGamma = _sidechain->findAtom("SG");
    AtomPtr hGamma = _sidechain->findAtom("HG");
    AtomPtr hBackBone = _sidechain->findAtom("HA");

    AtomPtr inherit = cAlpha;

    tieBetaCarbon(sGamma);

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

    AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
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

    tieBetaCarbon(cGamma1);

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

void Knotter::makeAlanine()
{
    MonomerPtr monomer = _sidechain->getMonomer();
    BackbonePtr backbone = monomer->getBackbone();
    std::string residue = monomer->getIdentifier();
    _sidechain->setCanRefine(true);

    AtomPtr spineAtom = backbone->betaCarbonTorsionAtom();
    AtomPtr cAlpha = _sidechain->findAtom("CA");
    AtomPtr cBeta = _sidechain->findAtom("CB");
    AtomPtr hBeta1 = _sidechain->findAtom("HB1");
    AtomPtr hBeta2 = _sidechain->findAtom("HB2");
    AtomPtr hBeta3 = _sidechain->findAtom("HB3");
    AtomPtr hBackbone = _sidechain->findAtom("HA");

    AtomPtr inherit = cAlpha;

    tieBetaCarbon(AtomPtr());

    BondPtr cb2hb1 = BondPtr(new Bond(cBeta, hBeta1));
    cb2hb1->activate(_sidechain, inherit);
    BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
    cb2hb2->activate(_sidechain, inherit);
    BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));
    cb2hb3->activate(_sidechain, inherit);
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


    BondPtr ca2cb = tieBetaCarbon(cGamma);
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

    cb2cg->setRefineBondAngle();
    
    BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
    cb2hb2->activate(_sidechain, inherit);
    BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));
    cb2hb3->activate(_sidechain, inherit);

    BondPtr cg2nd1 = BondPtr(new Bond(cGamma, nDelta1));
    cg2nd1->setFixed(true);
    cg2nd1->setTorsionAtoms(cBeta, cEpsilon1);
    cg2nd1->activate(_sidechain, inherit);

    BondPtr cg2cd2 = BondPtr(new Bond(cGamma, cDelta2));
    cg2cd2->setTorsionAtoms(cBeta, nEpsilon2);
    cg2cd2->setFixed(true);
    cg2cd2->activate(_sidechain, inherit);

    BondPtr nd12ce1 = BondPtr(new Bond(nDelta1, cEpsilon1));
    nd12ce1->setTorsionAtoms(cGamma, hEpsilon1);
    nd12ce1->setFixed(true);
    nd12ce1->activate(_sidechain, inherit);
    BondPtr nd12hd1 = BondPtr(new Bond(nDelta1, hDelta1));
    nd12hd1->activate(_sidechain, inherit);


    BondPtr ce12he1 = BondPtr(new Bond(cEpsilon1, hEpsilon1));
    ce12he1->activate(_sidechain, inherit);

    BondPtr cd22ne2 = BondPtr(new Bond(cDelta2, nEpsilon2));
    cd22ne2->setTorsionAtoms(cGamma, hEpsilon2);
    cd22ne2->setFixed(true);
    cd22ne2->activate(_sidechain, inherit);

    BondPtr cd22hd2 = BondPtr(new Bond(cDelta2, hDelta2));
    cd22hd2->activate(_sidechain, inherit);

    BondPtr ne22he2 = BondPtr(new Bond(nEpsilon2, hEpsilon2));
    ne22he2->activate(_sidechain, inherit);
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
    AtomPtr hBeta2 = _sidechain->findAtom("HB2");
    AtomPtr hBeta3 = _sidechain->findAtom("HB3");
    AtomPtr cGamma = _sidechain->findAtom("CG");
    AtomPtr cDelta1 = _sidechain->findAtom("CD1");
    AtomPtr hDelta1 = _sidechain->findAtom("HD1");
    AtomPtr cEpsilon1 = _sidechain->findAtom("CE1");
    AtomPtr hEpsilon1 = _sidechain->findAtom("HE1");
    AtomPtr cOmega = _sidechain->findAtom("CZ");
    AtomPtr oxygen = _sidechain->findAtom("OH");
    AtomPtr hydrogen = _sidechain->findAtom("HH");
    AtomPtr cEpsilon2 = _sidechain->findAtom("CE2");
    AtomPtr hEpsilon2 = _sidechain->findAtom("HE2");
    AtomPtr cDelta2 = _sidechain->findAtom("CD2");
    AtomPtr hDelta2 = _sidechain->findAtom("HD2");
    AtomPtr hBackbone = _sidechain->findAtom("HA");

    AtomPtr inherit = cAlpha;


    BondPtr ca2cb = tieBetaCarbon(cGamma);
    ca2cb->addExtraTorsionSample(cOmega, 0);
    ca2cb->addExtraTorsionSample(oxygen, 0);

    BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
    cb2cg->setTorsionAtoms(cAlpha, cDelta1);
    cb2cg->activate(_sidechain, inherit);
    cb2cg->addExtraTorsionSample(cEpsilon1, 0);
    cb2cg->addExtraTorsionSample(cEpsilon2, 0);

    ca2cb->setRefineBondAngle();
    cb2cg->setRefineBondAngle();

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
    ce22cz->setTorsionAtoms(cDelta2, oxygen);
    ce22cz->setFixed(true);
    ce22cz->activate(_sidechain, inherit);

    BondPtr ce22he2 = BondPtr(new Bond(cEpsilon2, hEpsilon2));
    ce22he2->activate(_sidechain, inherit);
    BondPtr cz2oh = BondPtr(new Bond(cOmega, oxygen));
    cz2oh->activate(_sidechain, inherit);
    BondPtr oh2hh = BondPtr(new Bond(oxygen, hydrogen));
    oh2hh->activate(_sidechain, inherit);
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


    BondPtr ca2cb = tieBetaCarbon(cGamma);
    if (!ca2cb) return;
    ca2cb->addExtraTorsionSample(cOmega, 0);

    BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
    cb2cg->setTorsionAtoms(cAlpha, cDelta1);
    cb2cg->activate(_sidechain, inherit);
    cb2cg->addExtraTorsionSample(cEpsilon1, 0);
    cb2cg->addExtraTorsionSample(cEpsilon2, 0);

    ca2cb->setRefineBondAngle();
    cb2cg->setRefineBondAngle();

    BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
    cb2hb2->activate(_sidechain, inherit);
    BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));
    cb2hb3->activate(_sidechain, inherit);

    BondPtr cg2cd1 = BondPtr(new Bond(cGamma, cDelta1));
    cg2cd1->setTorsionAtoms(cBeta, cEpsilon1);
    cg2cd1->setFixed(true);
    cg2cd1->activate(_sidechain, inherit);
    Bond::setTorsion(&*cg2cd1, deg2rad(180));

    BondPtr cg2cd2 = BondPtr(new Bond(cGamma, cDelta2));
    cg2cd2->setTorsionAtoms(cBeta, cEpsilon2);
    cg2cd2->setFixed(true);
    cg2cd2->activate(_sidechain, inherit);
    Bond::setTorsion(&*cg2cd2, deg2rad(180));

    BondPtr cd22ce2 = BondPtr(new Bond(cDelta2, cEpsilon2));
    cd22ce2->setTorsionAtoms(cGamma, cOmega);
    cd22ce2->setFixed(true);
    cd22ce2->activate(_sidechain, inherit);
    Bond::setTorsion(&*cd22ce2, deg2rad(0));
    BondPtr cd22hd2 = BondPtr(new Bond(cDelta2, hDelta2));
    cd22hd2->activate(_sidechain, inherit);

    BondPtr cd12ce1 = BondPtr(new Bond(cDelta1, cEpsilon1));
    cd12ce1->setTorsionAtoms(cGamma, hEpsilon1);
    cd12ce1->setFixed(true);
    cd12ce1->activate(_sidechain, inherit);
    Bond::setTorsion(&*cd12ce1, deg2rad(0));
    BondPtr cd12hd1 = BondPtr(new Bond(cDelta1, hDelta1));
    cd12hd1->activate(_sidechain, inherit);
    BondPtr ce12he1 = BondPtr(new Bond(cEpsilon1, hEpsilon1));
    ce12he1->activate(_sidechain, inherit);

    BondPtr ce22cz = BondPtr(new Bond(cEpsilon2, cOmega));
    ce22cz->setTorsionAtoms(cDelta2, hOmega);
    ce22cz->setFixed(true);
    ce22cz->activate(_sidechain, inherit);
    Bond::setTorsion(&*ce22cz, deg2rad(180));

    BondPtr ce22he2 = BondPtr(new Bond(cEpsilon2, hEpsilon2));
    ce22he2->activate(_sidechain, inherit);
    BondPtr cz2hz = BondPtr(new Bond(cOmega, hOmega));
    cz2hz->activate(_sidechain, inherit);

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
    AtomPtr hBeta2 = _sidechain->findAtom("HB2");
    AtomPtr hBeta3 = _sidechain->findAtom("HB3");
    AtomPtr cGamma = _sidechain->findAtom("CG");
    AtomPtr cDelta1 = _sidechain->findAtom("CD1");
    AtomPtr hDelta1 = _sidechain->findAtom("HD1");
    AtomPtr cDelta2 = _sidechain->findAtom("CD2");
    AtomPtr hDelta2 = _sidechain->findAtom("HD2");
    AtomPtr nEpsilon1 = _sidechain->findAtom("NE1");
    AtomPtr hEpsilon1 = _sidechain->findAtom("HE1");
    AtomPtr cEpsilon2 = _sidechain->findAtom("CE2");
    AtomPtr hEpsilon2 = _sidechain->findAtom("HE2");
    AtomPtr cEpsilon3 = _sidechain->findAtom("CE3");
    AtomPtr hEpsilon3 = _sidechain->findAtom("HE3");
    AtomPtr cOmega2 = _sidechain->findAtom("CZ2");
    AtomPtr hOmega2 = _sidechain->findAtom("HZ2");
    AtomPtr cOmega3 = _sidechain->findAtom("CZ3");
    AtomPtr hOmega3 = _sidechain->findAtom("HZ3");
    AtomPtr ch2 = _sidechain->findAtom("CH2");
    AtomPtr hh2 = _sidechain->findAtom("HH2");
    AtomPtr hBackbone = _sidechain->findAtom("HA");

    AtomPtr inherit = cAlpha;

    BondPtr ca2cb = tieBetaCarbon(cGamma);

    BondPtr cb2cg = BondPtr(new Bond(cBeta, cGamma));
    cb2cg->setTorsionAtoms(cAlpha, cDelta1);
    cb2cg->activate(_sidechain, inherit);
    cb2cg->addExtraTorsionSample(cDelta1, 0);
    cb2cg->addExtraTorsionSample(cDelta2, 0);
    cb2cg->addExtraTorsionSample(nEpsilon1, 0);
    cb2cg->addExtraTorsionSample(cEpsilon2, 0);
    cb2cg->addExtraTorsionSample(cEpsilon3, 0);
    cb2cg->addExtraTorsionSample(cOmega2, 0);
    cb2cg->addExtraTorsionSample(cOmega3, 0);
    cb2cg->addExtraTorsionSample(ch2, 0);

    ca2cb->setRefineBondAngle();
    cb2cg->setRefineBondAngle();

    BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
    cb2hb2->activate(_sidechain, inherit);
    BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));
    cb2hb3->activate(_sidechain, inherit);

    BondPtr cg2cd1 = BondPtr(new Bond(cGamma, cDelta1));
    cg2cd1->setTorsionAtoms(cBeta, nEpsilon1);
    cg2cd1->setFixed(true);
    cg2cd1->activate(_sidechain, inherit);

    BondPtr cg2cd2 = BondPtr(new Bond(cGamma, cDelta2));
    cg2cd2->setTorsionAtoms(cBeta, cEpsilon3);
    cg2cd2->setFixed(true);
    cg2cd2->activate(_sidechain, inherit);

    BondPtr cd22ce2 = BondPtr(new Bond(cDelta2, cEpsilon3));
    cd22ce2->setTorsionAtoms(cGamma, cOmega3);
    cd22ce2->setFixed(true);
    cd22ce2->activate(_sidechain, inherit);
//    BondPtr cd22hd2 = BondPtr(new Bond(cDelta2, hDelta2));
//    cd22hd2->activate(_sidechain, inherit);

    BondPtr cd12ce1 = BondPtr(new Bond(cDelta1, nEpsilon1));
    cd12ce1->setTorsionAtoms(cGamma, cEpsilon2);
    cd12ce1->setFixed(true);
    cd12ce1->activate(_sidechain, inherit);

    BondPtr cd12hd1 = BondPtr(new Bond(cDelta1, hDelta1));
    cd12hd1->activate(_sidechain, inherit);

    BondPtr ne12ce2 = BondPtr(new Bond(nEpsilon1, cEpsilon2));
    ne12ce2->setTorsionAtoms(cDelta1, cOmega2);
    ne12ce2->activate(_sidechain, inherit);

    BondPtr ce22cz = BondPtr(new Bond(cEpsilon2, cOmega2));
    ce22cz->setTorsionAtoms(cDelta2, hOmega2);
    ce22cz->setFixed(true);
    ce22cz->activate(_sidechain, inherit);

    BondPtr ce22he2 = BondPtr(new Bond(cEpsilon3, cOmega3));
    ce22he2->setTorsionAtoms(cDelta2, ch2);
    ce22he2->activate(_sidechain, inherit);
    BondPtr cz2hz = BondPtr(new Bond(cOmega3, hOmega3));
    cz2hz->activate(_sidechain, inherit);

    BondPtr ce32he3 = BondPtr(new Bond(cOmega3, ch2));
    ce32he3->activate(_sidechain, inherit);
    BondPtr ch2hh2 = BondPtr(new Bond(ch2, hh2));
    ch2hh2->activate(_sidechain, inherit);

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
    AtomPtr hBeta = _sidechain->findAtom("HB");
    AtomPtr cGamma1 = _sidechain->findAtom("CG1");
    AtomPtr hGamma12 = _sidechain->findAtom("HG12");
    AtomPtr hGamma13 = _sidechain->findAtom("HG13");
    AtomPtr cGamma2 = _sidechain->findAtom("CG2");
    AtomPtr cDelta1 = _sidechain->findAtom("CD1");
    AtomPtr hGamma21 = _sidechain->findAtom("HG21");
    AtomPtr hGamma22 = _sidechain->findAtom("HG22");
    AtomPtr hGamma23 = _sidechain->findAtom("HG23");
    AtomPtr hDelta11 = _sidechain->findAtom("HD11");
    AtomPtr hDelta12 = _sidechain->findAtom("HD12");
    AtomPtr hDelta13 = _sidechain->findAtom("HD13");
    AtomPtr hBackbone = _sidechain->findAtom("HA");
    AtomPtr inherit = cAlpha;

    BondPtr ca2cb = tieBetaCarbon(cGamma1);

    BondPtr cb2cg1 = BondPtr(new Bond(cBeta, cGamma1));
    cb2cg1->setTorsionAtoms(cAlpha, cDelta1);
    cb2cg1->activate(_sidechain, inherit);
    BondPtr cb2hb = BondPtr(new Bond(cBeta, hBeta));
    cb2hb->activate(_sidechain, inherit);
    BondPtr cb2cg2 = BondPtr(new Bond(cBeta, cGamma2));
    cb2cg2->activate(_sidechain, inherit);

    BondPtr cg1cd11 = BondPtr(new Bond(cGamma1, cDelta1));
    cg1cd11->activate(_sidechain, inherit);
    BondPtr cg1hg12 = BondPtr(new Bond(cGamma1, hGamma12));
    cg1hg12->activate(_sidechain, inherit);
    BondPtr cg1hg13 = BondPtr(new Bond(cGamma1, hGamma13));
    cg1hg13->activate(_sidechain, inherit);

    BondPtr cg2hg21 = BondPtr(new Bond(cGamma2, hGamma21));
    cg2hg21->activate(_sidechain, inherit);
    BondPtr cg2hg22 = BondPtr(new Bond(cGamma2, hGamma22));
    cg2hg22->activate(_sidechain, inherit);
    BondPtr cg2hg23 = BondPtr(new Bond(cGamma2, hGamma23));
    cg2hg23->activate(_sidechain, inherit);

    BondPtr cd1hg21 = BondPtr(new Bond(cDelta1, hDelta11));
    cd1hg21->activate(_sidechain, inherit);
    BondPtr cd1hg22 = BondPtr(new Bond(cDelta1, hDelta12));
    cd1hg22->activate(_sidechain, inherit);
    BondPtr cd1hg23 = BondPtr(new Bond(cDelta1, hDelta13));
    cd1hg23->activate(_sidechain, inherit);

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
    AtomPtr hBeta2 = _sidechain->findAtom("HB2");
    AtomPtr hBeta3 = _sidechain->findAtom("HB3");
    AtomPtr cGamma = _sidechain->findAtom("CG");
    AtomPtr hGamma = _sidechain->findAtom("HG");
    AtomPtr cDelta1 = _sidechain->findAtom("CD1");
    AtomPtr hDelta11 = _sidechain->findAtom("HD11");
    AtomPtr hDelta12 = _sidechain->findAtom("HD12");
    AtomPtr hDelta13 = _sidechain->findAtom("HD13");
    AtomPtr cDelta2 = _sidechain->findAtom("CD2");
    AtomPtr hDelta21 = _sidechain->findAtom("HD21");
    AtomPtr hDelta22 = _sidechain->findAtom("HD22");
    AtomPtr hDelta23 = _sidechain->findAtom("HD23");
    AtomPtr hBackbone = _sidechain->findAtom("HA");

    AtomPtr inherit = cAlpha;

    BondPtr ca2cb = tieBetaCarbon(cGamma);

    BondPtr cb2cg1 = BondPtr(new Bond(cBeta, cGamma));
    cb2cg1->setTorsionAtoms(cAlpha, cDelta1);
    cb2cg1->activate(_sidechain, inherit);
    BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
    cb2hb2->activate(_sidechain, inherit);
    BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));
    cb2hb3->activate(_sidechain, inherit);

    BondPtr cg1cd1 = BondPtr(new Bond(cGamma, cDelta1));
    cg1cd1->activate(_sidechain, inherit);
    BondPtr cg1cd2 = BondPtr(new Bond(cGamma, cDelta2));
    cg1cd2->activate(_sidechain, inherit);
    BondPtr cg1hg = BondPtr(new Bond(cGamma, hGamma));
    cg1hg->activate(_sidechain, inherit);

    BondPtr cd1hd11 = BondPtr(new Bond(cDelta1, hDelta11));
    cd1hd11->activate(_sidechain, inherit);
    BondPtr cd1hd12 = BondPtr(new Bond(cDelta1, hDelta12));
    cd1hd12->activate(_sidechain, inherit);
    BondPtr cd1hd13 = BondPtr(new Bond(cDelta1, hDelta13));
    cd1hd13->activate(_sidechain, inherit);

    BondPtr cd2hd21 = BondPtr(new Bond(cDelta2, hDelta21));
    cd2hd21->activate(_sidechain, inherit);
    BondPtr cd2hd22 = BondPtr(new Bond(cDelta2, hDelta22));
    cd2hd22->activate(_sidechain, inherit);
    BondPtr cd2hd23 = BondPtr(new Bond(cDelta2, hDelta23));
    cd2hd23->activate(_sidechain, inherit);
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
    AtomPtr hBeta2 = _sidechain->findAtom("HB2");
    AtomPtr hBeta3 = _sidechain->findAtom("HB3");
    AtomPtr cGamma = _sidechain->findAtom("CG");
    AtomPtr oDelta1 = _sidechain->findAtom("OD1");
    AtomPtr oDelta2 = _sidechain->findAtom("OD2");
    AtomPtr hBackbone = _sidechain->findAtom("HA");

    AtomPtr inherit = cAlpha;

    BondPtr ca2cb = tieBetaCarbon(cGamma);

    BondPtr cb2cg1 = BondPtr(new Bond(cBeta, cGamma));
    cb2cg1->setTorsionAtoms(cAlpha, oDelta1);
    cb2cg1->activate(_sidechain, inherit);

    cb2cg1->setRefineBondAngle();

    BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
    cb2hb2->activate(_sidechain, inherit);
    BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));
    cb2hb3->activate(_sidechain, inherit);

    BondPtr cg1cd1 = BondPtr(new Bond(cGamma, oDelta1));
    cg1cd1->activate(_sidechain, inherit);
    BondPtr cg1cd2 = BondPtr(new Bond(cGamma, oDelta2));
    cg1cd2->activate(_sidechain, inherit);    
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
    AtomPtr hBeta2 = _sidechain->findAtom("HB2");
    AtomPtr hBeta3 = _sidechain->findAtom("HB3");
    AtomPtr cGamma = _sidechain->findAtom("CG");
    AtomPtr oDelta1 = _sidechain->findAtom("OD1");
    AtomPtr nDelta2 = _sidechain->findAtom("ND2");
    AtomPtr hDelta21 = _sidechain->findAtom("HD21");
    AtomPtr hDelta22 = _sidechain->findAtom("HD22");
    AtomPtr hBackbone = _sidechain->findAtom("HA");

    AtomPtr inherit = cAlpha;

    BondPtr ca2cb = tieBetaCarbon(cGamma);

    BondPtr cb2cg1 = BondPtr(new Bond(cBeta, cGamma));
    cb2cg1->setTorsionAtoms(cAlpha, oDelta1);
    cb2cg1->activate(_sidechain, inherit);
    BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
    cb2hb2->activate(_sidechain, inherit);
    BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));
    cb2hb3->activate(_sidechain, inherit);

    BondPtr cg1nd2 = BondPtr(new Bond(cGamma, nDelta2));
    cg1nd2->activate(_sidechain, inherit);
    BondPtr cg1od1 = BondPtr(new Bond(cGamma, oDelta1));
    cg1od1->activate(_sidechain, inherit);
    BondPtr nd22hd21 = BondPtr(new Bond(nDelta2, hDelta21));
    nd22hd21->activate(_sidechain, inherit);
    BondPtr nd22hd22 = BondPtr(new Bond(nDelta2, hDelta22));
    nd22hd22->activate(_sidechain, inherit);
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
    AtomPtr hBeta2 = _sidechain->findAtom("HB2");
    AtomPtr hBeta3 = _sidechain->findAtom("HB3");
    AtomPtr cGamma = _sidechain->findAtom("CG");
    AtomPtr hGamma2 = _sidechain->findAtom("HG2");
    AtomPtr hGamma3 = _sidechain->findAtom("HG3");
    AtomPtr cDelta = _sidechain->findAtom("CD");
    AtomPtr oEpsilon1 = _sidechain->findAtom("OE1");
    AtomPtr nEpsilon2 = _sidechain->findAtom("NE2");
    AtomPtr hEpsilon21 = _sidechain->findAtom("HE21");
    AtomPtr hEpsilon22 = _sidechain->findAtom("HE22");
    AtomPtr hBackbone = _sidechain->findAtom("HA");

    AtomPtr inherit = cAlpha;

    BondPtr ca2cb = tieBetaCarbon(cGamma);

    BondPtr cb2cg1 = BondPtr(new Bond(cBeta, cGamma));
    cb2cg1->setTorsionAtoms(cAlpha, cDelta);
    cb2cg1->activate(_sidechain, inherit);
    BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
    cb2hb2->activate(_sidechain, inherit);
    BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));
    cb2hb3->activate(_sidechain, inherit);

    BondPtr cg2cd1 = BondPtr(new Bond(cGamma, cDelta));
    cg2cd1->setTorsionAtoms(cAlpha, oEpsilon1);
    cg2cd1->activate(_sidechain, inherit);
    BondPtr cg2hg2 = BondPtr(new Bond(cGamma, hGamma2));
    cg2hg2->activate(_sidechain, inherit);
    BondPtr cg2hg3 = BondPtr(new Bond(cGamma, hGamma3));
    cg2hg3->activate(_sidechain, inherit);

    BondPtr cd2ce1 = BondPtr(new Bond(cDelta, oEpsilon1));
    cd2ce1->activate(_sidechain, inherit);
    BondPtr cd2ce2 = BondPtr(new Bond(cDelta, nEpsilon2));
    cd2ce2->activate(_sidechain, inherit);
    BondPtr nd22hd21 = BondPtr(new Bond(nEpsilon2, hEpsilon21));
    nd22hd21->activate(_sidechain, inherit);
    BondPtr nd22hd22 = BondPtr(new Bond(nEpsilon2, hEpsilon22));
    nd22hd22->activate(_sidechain, inherit);
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
    AtomPtr hBeta2 = _sidechain->findAtom("HB2");
    AtomPtr hBeta3 = _sidechain->findAtom("HB3");
    AtomPtr cGamma = _sidechain->findAtom("CG");
    AtomPtr hGamma2 = _sidechain->findAtom("HG2");
    AtomPtr hGamma3 = _sidechain->findAtom("HG3");
    AtomPtr cDelta = _sidechain->findAtom("CD");
    AtomPtr oEpsilon1 = _sidechain->findAtom("OE1");
    AtomPtr oEpsilon2 = _sidechain->findAtom("OE2");
    AtomPtr hBackbone = _sidechain->findAtom("HA");

    AtomPtr inherit = cAlpha;

    BondPtr ca2cb = tieBetaCarbon(cGamma);

    BondPtr cb2cg1 = BondPtr(new Bond(cBeta, cGamma));
    cb2cg1->setTorsionAtoms(cAlpha, cDelta);
    cb2cg1->activate(_sidechain, inherit);
    BondPtr cb2hb2 = BondPtr(new Bond(cBeta, hBeta2));
    cb2hb2->activate(_sidechain, inherit);
    BondPtr cb2hb3 = BondPtr(new Bond(cBeta, hBeta3));
    cb2hb3->activate(_sidechain, inherit);

    BondPtr cg2cd1 = BondPtr(new Bond(cGamma, cDelta));
    cg2cd1->setTorsionAtoms(cAlpha, oEpsilon1);
    cg2cd1->activate(_sidechain, inherit);
    BondPtr cg2hg2 = BondPtr(new Bond(cGamma, hGamma2));
    cg2hg2->activate(_sidechain, inherit);
    BondPtr cg2hg3 = BondPtr(new Bond(cGamma, hGamma3));
    cg2hg3->activate(_sidechain, inherit);

    BondPtr cd2ce1 = BondPtr(new Bond(cDelta, oEpsilon1));
    cd2ce1->activate(_sidechain, inherit);
    BondPtr cd2ce2 = BondPtr(new Bond(cDelta, oEpsilon2));
    cd2ce2->activate(_sidechain, inherit);
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
    AtomPtr hBeta = _sidechain->findAtom("HB");
    AtomPtr oGamma1 = _sidechain->findAtom("OG1");
    AtomPtr hGamma11 = _sidechain->findAtom("HG1");
    AtomPtr cGamma2 = _sidechain->findAtom("CG2");
    AtomPtr hGamma21 = _sidechain->findAtom("HG21");
    AtomPtr hGamma22 = _sidechain->findAtom("HG22");
    AtomPtr hGamma23 = _sidechain->findAtom("HG23");
    AtomPtr hBackbone = _sidechain->findAtom("HA");

    AtomPtr inherit = cAlpha;

    BondPtr ca2cb = tieBetaCarbon(oGamma1);

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
