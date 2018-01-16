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
	_strategy->setCycles(40);
//	_strategy->setVerbose(true);

	_strategy->setJobName("flex_region");
	_strategy->setSilent(false);
}

void FlexRegion::sample()
{
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

	_strategy->addParameter(&*_bonds[i], Bond::getTorsionBlur, Bond::setTorsionBlur, 0.001, 0.0001, "kick");
	_strategy->addParameter(&*_bonds[i], Bond::getDampening, Bond::setDampening, 0.005, 0.0001, "dampen");
//	_strategy->addParameter(&*_bonds[i], Bond::getMagicPhi, Bond::setMagicPhi, deg2rad(10), deg2rad(1), "phi");
//	_strategy->addParameter(&*_bonds[i], Bond::getMagicPsi, Bond::setMagicPsi, deg2rad(10), deg2rad(1), "psi");
}

double FlexRegion::getScore()
{
	if (!bondCount()) return 0;

	double score = 0;
	double count = 0;

	_bonds[0]->propagateChange(bondCount());
	MoleculePtr molecule = _bonds[0]->getMolecule();
	PolymerPtr polymer = ToPolymerPtr(molecule);

	ModelPtr model = polymer->getAnchorModel();
	mat3x3 baseTensor = model->getRealSpaceTensor();
	mat3x3 origBaseTensor = make_mat3x3();
	mat3x3 inverseBase = mat3x3_inverse(baseTensor);

	if (model->isAbsolute())
	{
		origBaseTensor = ToAbsolutePtr(model)->getAtom()->getTensor();
	}
	else
	{
		shout_at_helen("Anchor is not an Absolute!");
	}

	mat3x3 baseDiff = mat3x3_mult_mat3x3(origBaseTensor, inverseBase);

	for (int i = 0; i < bondCount(); i++)
	{
		BondPtr checkBond = _bonds[i];
		
		AtomPtr minor = checkBond->getMinor();
		AtomPtr reference = ToAbsolutePtr(model)->getAtom();

		/* Dealing with our Vagabondage distributions */

		mat3x3 actualRefTensor = model->getRealSpaceTensor();
		mat3x3 actualRefRemainder = mat3x3_mult_mat3x3(baseDiff, actualRefTensor);
		mat3x3 myActualTensor = checkBond->getRealSpaceTensor();
		mat3x3 myActualRemainder = mat3x3_mult_mat3x3(baseDiff, myActualTensor);
		mat3x3 myActualInverse = mat3x3_inverse(myActualRemainder);
		mat3x3 actualDifference = mat3x3_mult_mat3x3(actualRefRemainder, myActualInverse);
	//	mat3x3 actualInverse = mat3x3_inverse(actualDifference);

	//	mat3x3 mult = mat3x3_mult_mat3x3(actualInverse, idealDifference);
	//	mat3x3 diff_inverse = mat3x3_inverse(bond_diff);

	//	mat3x3 mult = mat3x3_mult_mat3x3(diff_inverse, difference);

		// Mult should be identity if it is perfect...

		double diff = mat3x3_diff_from_identity(actualDifference, 1);
		score += diff;
		count++;
	}

	return score / count;
}
