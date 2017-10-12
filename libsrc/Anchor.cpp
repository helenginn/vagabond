//
//  Anchor.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/09/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Anchor.h"
#include "Bond.h"
#include "Shouter.h"

Anchor::Anchor(BondPtr inheritToNTerminus, BondPtr inheritToCTerminus)
{
	/* When going from towards-C-tie (>) to towards-N-tie (<):
	 *                *
	 *    --N--CA--C--N--CA--C--N--CA--
	 *              ^   ^
	 *    downstream     parent
	 *
	 * This class must describe *N in a way to prepare the preceding C
	 * i.e. start = position of N, old_start is position of CA,
	 * basis is CA--N bond reversed, torsion is reversed.
	 */

	_trappedToNTerminus = std::make_shared<Bond>(*inheritToNTerminus);
	_trappedToCTerminus = std::make_shared<Bond>(*inheritToCTerminus);

	// Either we have bonds C->N and N->CA (to C terminus)
	// or we have CA->N and N->C (to N terminus)
	// We must flip accordingly.
	_flipNTerminus = (_trappedToCTerminus->getMinor()->getAtomName() == "N");

	BondPtr badBond = _trappedToCTerminus;
	BondPtr goodBond = _trappedToNTerminus;



	if (_flipNTerminus)
	{
		badBond = _trappedToNTerminus;
		goodBond = _trappedToCTerminus;
		_flipNTerminus = true;
	}

	_anchorStart = goodBond->getStaticPosition();
	_newHeavyStoredGroup = ToBondPtr(badBond->downstreamAtom(0, 0)->getModel())->getBondGroup(0);

	// we need to flip the bond to the N-terminus
	badBond->reverse(BondPtr());

	if (badBond->downstreamAtomGroupCount() != 1)
	{
		shout_at_helen("Helen cannot anchor a bond which has two downstream groups");
	}

	badBond->getBondGroup(0)->atoms.clear();

	/* Add all the upstream bonds as downstream bonds during semi-manual reversal */
	badBond->addDownstreamAtom(goodBond->getMajor(), 0);

	for (int j = goodBond->downstreamAtomCount(0) - 1; j > 0; j--)
	{
		AtomPtr atom = goodBond->downstreamAtom(0, j);
		badBond->addDownstreamAtom(atom, 0);
	}

	if (_trappedToNTerminus->getMinor() != _trappedToCTerminus->getMinor())
	{
		std::cout << "Fuck, it's " << _trappedToNTerminus->shortDesc()
		<< " and " << _trappedToCTerminus->shortDesc() << std::endl;
	}
}

void Anchor::activate()
{
	// This tells the bonds not to respond to propagation of changes anymore.
	_trappedToCTerminus->setAnchored();
	_trappedToNTerminus->setAnchored();

	BondPtr badBond = _trappedToCTerminus;
	BondPtr goodBond = _trappedToNTerminus;

	if (_flipNTerminus) // had to flip the left hand side to face N terminus
	{
		badBond = _trappedToNTerminus; // bad bond is the flipped one.
		goodBond = _trappedToCTerminus;
	}
	BondGroup *badGroup = badBond->getBondGroup(0);
	BondGroup *goodGroup = goodBond->getBondGroup(0);

	mat3x3 reverseDir = make_mat3x3();
	vec3 x = make_vec3(1, 0, 0);
	vec3 z = make_vec3(0, 0, -1);
	reverseDir = mat3x3_rhbasis(x, z);

	BondSample *staticSample = &badGroup->staticSample[0];
	staticSample->start = _anchorStart;

	for (int j = 0; j < badGroup->storedSamples.size(); j++)
	{
		BondSample *preSample = &_newHeavyStoredGroup->storedSamples[j];
		BondSample *sample = &badGroup->storedSamples[j];
		BondSample *goodSample = &goodGroup->storedSamples[j];

		vec3 lPos = goodSample->old_start;
		vec3 miPos = goodSample->start;
		vec3 maPos = sample->start;
		vec3 hPos = preSample->start;

		sample->basis = badBond->makeTorsionBasis(hPos, maPos, miPos, lPos,
												  &sample->torsion);
		sample->old_start = sample->start;
		sample->start = goodSample->start;

	}

	AtomPtr nitrogen = _trappedToNTerminus->getMinor();
	nitrogen->setModel(shared_from_this());
}

BondPtr Anchor::getAppropriateBond()
{
	if (_trappedToCTerminus->getMajor() == _callingBond->getMinor())
	{
		return _trappedToNTerminus;
	}

	return _trappedToCTerminus;
}

BondPtr Anchor::sanitiseBond(Bond *myself, BondPtr model)
{
	if (model->isAnchor())
	{
		AnchorPtr anchor = ToAnchorPtr(ToBondPtr(model));
		anchor->setCallingBond(myself);
		model = anchor->getAppropriateBond();
	}

	return model;
}
