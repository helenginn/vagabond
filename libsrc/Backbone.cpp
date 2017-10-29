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
	return (bond->getMinor()->getAtomName() == "CA" ||
			bond->getMajor()->getAtomName() == "CA");
}

void Backbone::refine(CrystalPtr target, RefinementType rType)
{
	if (!isTied())
	{
		return;
	}

	std::cout << getMonomer()->getResCode() << std::flush;

	AtomGroup::refine(target, rType);
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
