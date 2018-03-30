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
#include "Sidechain.h"

bool Backbone::shouldRefineMagicAxis(BondPtr)
{
	return false;
}

void Backbone::refine(CrystalPtr target, RefinementType rType)
{
	if (!isTied())
	{
		return;
	}
	
	if (rType == RefinementSidechain)
	{
		return;
	}

	std::cout << getMonomer()->getResCode() << std::flush;

	if (!paramCount())
	{
		double range = 2.;

		if (_timesRefined >= 3)
		{
			range = 0.1;
		}

		switch (rType)
		{
			case RefinementModelPos:
			addParamType(ParamOptionTorsion, range);
			break;

			case RefinementFine:
			addParamType(ParamOptionTorsion, range);
			break;

			default:
			break;
		}
	}

	MonomerPtr monomer = getMonomer();
	SidechainPtr sidechain = monomer->getSidechain();
	if (sidechain && rType == RefinementFine)
	{
//		addIncludeForRefinement(sidechain);
	}

	AtomGroup::refine(target, rType);
	clearParams();
	clearIncludeForRefinements();
}

AtomPtr Backbone::betaCarbonTorsionAtom()
{
	/* What is the major atom of the bond describing CA? */

	AtomPtr ca = findAtom("CA");

	if (!ca)
	{
		return AtomPtr();
	}

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
	// This will become dodgy if we have multiple conformers on backbone
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
		for (int k = 0; k < bond->downstreamAtomGroupCount(); k++)
		{
			oldReversal->getBondGroup(k)->atoms.clear();

			if (currentReversal->isAnchor())
			{
				ToAnchorPtr(currentReversal)->setCallingBond(&*oldReversal);
				BondPtr bond = ToAnchorPtr(currentReversal)->getAppropriateBond(true);

				oldReversal->addDownstreamAtom(bond->getMajor(), k);

				for (int i = bond->downstreamAtomCount(k) - 1; i > 0; i--)
				{
					oldReversal->addDownstreamAtom(bond->downstreamAtom(k, i), 0);
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
						oldReversal->addDownstreamAtom(atom, k);
					}
				}
			}
		}
	}
}

void Backbone::addProperties()
{
	addIntProperty("res_num", &_resNum);
	AtomGroup::addProperties();
}
