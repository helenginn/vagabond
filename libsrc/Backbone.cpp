//
//  Backbone.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Backbone.h"
#include "Monomer.h"
#include "Bond.h"
#include "Atom.h"
#include "Absolute.h"
#include "FileReader.h"
#include "Anchor.h"

bool Backbone::shouldRefineMagicAxis(BondPtr bond)
{
//	if (_refinedMagicAxisCount > 2) return false;

	return (bond->getMinor()->getAtomName() == "CA" ||
			bond->getMajor()->getAtomName() == "CA");
}

void Backbone::refine(CrystalPtr target, RefinementType rType)
{
	int resNum = getMonomer()->getResidueNum();

	if (!isTied())
	{
		return;
	}

	std::cout << getMonomer()->getResCode() << std::flush;

	for (int i = 0; i < atomCount(); i++)
	{
		AtomPtr myAtom = atom(i);
		ModelPtr model = myAtom->getModel();

		if (!model->isBond())
		{
			continue;
		}

		BondPtr bond = ToBondPtr(model);

		if (!bond->isRefinable())
		{
			continue;
		}

		int groups = bond->downstreamAtomGroupCount();

		for (int k = 0; k < groups; k++)
		{
			if (rType != RefinementModelRMSD)
			{
				break;
			}

			int magicStart = resNum;
			int magicEnd = resNum + FUTURE_RESIDUES;
			int magicCloseEnd = resNum + 2;

			if (!getMonomer()->isAfterAnchor())
			{
				magicEnd = resNum - FUTURE_RESIDUES;
				magicCloseEnd = resNum - 2;
			}

			if (shouldRefineMagicAxis(bond))
			{
				_refinedMagicAxisCount++;

				/* Find the vague direction */
				bond->calculateMagicAxis();

				/* Fine-tune the axis */
			/*	setupGrid();
				addMagicAngle(bond, deg2rad(360.0), deg2rad(60.0));
				addSampledBackbone(getPolymer(), magicStart, magicEnd);
				setJobName("magic_angle_" + bond->shortDesc());
				setSilent();
				setScoreType(ScoreTypeModelRMSDZero);
				sample();

				setupNelderMead();
				setCycles(6);
				addMagicAngle(bond, deg2rad(15.0), deg2rad(2.0));
				addSampledBackbone(getPolymer(), magicStart, magicEnd);
				setJobName("magic_angle_" + bond->shortDesc());
				setSilent();
				setScoreType(ScoreTypeModelRMSDZero);
				sample();*/

			//	double angle = Bond::getMagicAngle(&*bond);
			//	std::cout << std::endl << "angle\t" << angle << std::endl;

			}

			/* Refine model position */

			setupNelderMead();

			if (getMonomer()->isAfterAnchor())
			{
				addRamachandranAngles(getPolymer(), resNum, magicCloseEnd);
				addSampledBackbone(getPolymer(), resNum, magicCloseEnd);
			}
			else
			{
				addRamachandranAngles(getPolymer(), resNum, magicCloseEnd);
				addSampledBackbone(getPolymer(), resNum, magicCloseEnd);
			}

			setScoreType(ScoreTypeModelPos);
			setSilent();
			setJobName("model_pos_" +  bond->shortDesc());
			sample();
		}

		if (rType != RefinementFineBlur)
		{
			continue;
		}

		for (int k = 0; k < groups; k++)
		{
			setupNelderMead();
			setupTorsionSet(bond, k, 5, resNum, 0.05, 0.02);
			setCrystal(target);
			sample();

			setupNelderMead();
			setupTorsionSet(bond, k, 4, resNum, 0.3, 0.02);
			setCrystal(target);
			sample();
		}

		bond->setActiveGroup(0);
	}
}

AtomPtr Backbone::betaCarbonTorsionAtom()
{
	/* What is the major atom of the bond describing CA? */

	AtomPtr ca = findAtom("CA");
	ModelPtr model = ca->getModel();

	if (model->getClassName() == "Bond")
	{
		BondPtr bond = ToBondPtr(model);
		return bond->getMajor();
	}

	return AtomPtr();
}

void Backbone::setAnchor()
{
	// the backbone N should become an anchor point
	AtomPtr nitrogen = findAtom("N");
	ModelPtr model = nitrogen->getModel();
	BondPtr bond = ToBondPtr(model);
	AtomPtr downstreamAtom = bond->downstreamAtom(0, 0);
	ModelPtr nextModel = downstreamAtom->getModel();
	BondPtr nextBond = ToBondPtr(nextModel);

	AnchorPtr anchor = AnchorPtr(new Anchor(bond, nextBond));
	ModelPtr modelAnchor = ToModelPtr(ToBondPtr(anchor));

	ModelPtr currentReversal = model;
	BondPtr oldReversal = BondPtr();

	while (currentReversal->getClassName() == "Bond")
	{
		BondPtr reverseBond = ToBondPtr(currentReversal);
		currentReversal = reverseBond->reverse(oldReversal);
		oldReversal = reverseBond;
	}

	std::cout << "Reanchoring on atom " << nitrogen->shortDesc() << std::endl;

	ToAnchorPtr(modelAnchor)->activate();

	if (currentReversal->getClassName() == "Absolute" ||
		currentReversal->isAnchor())
	{
		oldReversal->getBondGroup(0)->atoms.clear();

		if (currentReversal->isAnchor())
		{
			ToAnchorPtr(currentReversal)->setCallingBond(&*oldReversal);
			BondPtr bond = ToAnchorPtr(currentReversal)->getAppropriateBond(true);

			oldReversal->addDownstreamAtom(bond->getMajor(), 0);

			for (int i = bond->downstreamAtomCount(0) - 1; i > 0; i--)
			{
				oldReversal->addDownstreamAtom(bond->downstreamAtom(0, i), 0);
			}
		}
		else
		{
			AbsolutePtr absolute = ToAbsolutePtr(currentReversal);
			for (int i = 0; i < absolute->nextAtomCount(); i++)
			{
				AtomPtr atom = absolute->getNextAtom(i);
				if (atom != oldReversal->getMajor())
				{
					oldReversal->addDownstreamAtom(atom, 0);
				}
			}
		}
	}
}
