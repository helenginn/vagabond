//
//  FlexRegion.cpp
//  vagabond
//
//  Created by Helen Ginn on 04/01/2018.
//  Copyright © 2018 Strubi. All rights reserved.
//

#include "FlexRegion.h"
#include "Bond.h"
#include "Absolute.h"
#include <algorithm>
#include "Shouter.h"
#include "Polymer.h"
#include "RefinementNelderMead.h"

FlexRegion::FlexRegion()
{
    
}

void FlexRegion::addBond(BondPtr bond, int prevBondCount)
{
    _bonds.clear();
    _bonds.push_back(bond);
    BondPtr lastBond = bond;

    while (prevBondCount > 0)
    {
        ModelPtr model = lastBond->getParentModel();
        if (model && model->isBond())
        {
            lastBond = ToBondPtr(model);
            _bonds.push_back(lastBond);
            prevBondCount--;
        }
        else
        {
            break;
        }
    }

    std::reverse(_bonds.begin(), _bonds.end());
}

void FlexRegion::setup()
{
    NelderMeadPtr mead = NelderMeadPtr(new NelderMead());
    _strategy = boost::static_pointer_cast<RefinementStrategy>(mead);
    _strategy->setEvaluationFunction(FlexRegion::score, this);
    _strategy->setCycles(100);

    _strategy->setJobName("flex_region");
    _strategy->setSilent(false);
}

void FlexRegion::sample()
{
    if (!bondCount())
    {
        _bonds.clear();
        _strategy = RefinementStrategyPtr();
        return;
    }

    _strategy->setJobName("flex_region_" + _bonds[0]->shortDesc());
    _strategy->refine();
}

void FlexRegion::addSingleBondParameters()
{
    for (int i = 0; i < bondCount() - 2; i++)
    {
        addSingleBondParameter(i);
    }
}

void FlexRegion::addSingleBondParameter(int i)
{
    if (!_strategy)
    {
        shout_at_helen("Trying to add single bond param when strategy not ready.");
    }

    _strategy->addParameter(&*_bonds[i], Bond::getTorsionBlur, Bond::setTorsionBlur, 0.002, 0.0001, "kick");
    _strategy->addParameter(&*_bonds[i], Bond::getDampening, Bond::setDampening, 0.01, 0.0001, "dampen");
    _strategy->addParameter(&*_bonds[i], Bond::getMagicPhi, Bond::setMagicPhi, deg2rad(10), deg2rad(1), "phi");
    _strategy->addParameter(&*_bonds[i], Bond::getMagicPsi, Bond::setMagicPsi, deg2rad(10), deg2rad(1), "psi");
}

double FlexRegion::getScore()
{
    if (!bondCount()) return 0;

    double score = 0;
    double count = 0;

    _bonds[0]->propagateChange(bondCount());
    MoleculePtr molecule = _bonds[0]->getMolecule();
    PolymerPtr polymer = ToPolymerPtr(molecule);
    ModelPtr anchor = polymer->getAnchorModel();

    /* Current tensor for vagabond model */
    mat3x3 baseTensor = anchor->getRealSpaceTensor();

    /* Original tensor from the PDB */
    mat3x3 origBaseTensor = make_mat3x3();

    if (anchor->isAbsolute())
    {
        origBaseTensor = ToAbsolutePtr(anchor)->getAtom()->getTensor();
    }
    else
    {
        shout_at_helen("Anchor is not an Absolute!");
    }

    mat3x3 constantDiff = mat3x3_subtract_mat3x3(origBaseTensor, baseTensor);

    /* Calculate inflation from 'current' to 'original' */

    for (int i = 0; i < bondCount(); i++)
    {
        BondPtr checkBond = _bonds[i];
        AtomPtr checkAtom = checkBond->getMinor();

        /* Take the current tensor and inflate it to match anchor inflation */
        mat3x3 currentTensor = checkBond->getRealSpaceTensor();

        /* Take the original tensor which should be somewhat inflated already */
        mat3x3 origTensor = checkAtom->getTensor();

        /* This bit is actually now active */
        mat3x3 removeConstant = mat3x3_subtract_mat3x3(origTensor,
                constantDiff);

        mat3x3 finalDiff = mat3x3_subtract_mat3x3(currentTensor, removeConstant);

        /* Find how close this is to the identity matrix */
        double diff = mat3x3_abs_sum_all(finalDiff);
        score += diff;
        count++;
    }

    return score / count;
    /* Version comparing bond by bond */

    for (int i = 0; i < bondCount() - 1; i++)
    {
        for (int j = 0; j < i; j++)
        {
            BondPtr lastBond = _bonds[i];
            AtomPtr lastAtom = lastBond->getMinor();

            BondPtr checkBond = _bonds[j];
            AtomPtr checkAtom = checkBond->getMinor();

            /* Take the current tensor and inflate it to match anchor inflation */
            mat3x3 currentTensor = checkBond->getRealSpaceTensor();

            /* Take the original tensor which should be somewhat inflated already */
            mat3x3 origTensor = checkAtom->getTensor();

            mat3x3 lastCurrent = lastBond->getRealSpaceTensor();
            mat3x3 lastOrig = lastAtom->getTensor();

            mat3x3 aimDiff = mat3x3_subtract_mat3x3(origTensor, lastOrig);
            mat3x3 ourDiff = mat3x3_subtract_mat3x3(currentTensor, lastCurrent);

            mat3x3 finalDiff = mat3x3_subtract_mat3x3(ourDiff, aimDiff);

            /* Find how close this is to the identity matrix */
            double diff = mat3x3_abs_sum_all(finalDiff);
            score += diff;
            count++;
        }
    }

    return score / count;
}
