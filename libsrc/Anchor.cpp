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

	_callingBond = NULL;
	_trappedToNTerminus = std::make_shared<Bond>(*inheritToNTerminus);
	_trappedToCTerminus = std::make_shared<Bond>(*inheritToCTerminus);

	// Either we have bonds C->N and N->CA (to C terminus)
	// or we have CA->N and N->C (to N terminus)
	// We must flip accordingly. First we determine the state.
	_flipNTerminus = (_trappedToCTerminus->getMinor()->getAtomName() == "N");

	// "Bad bond" is the one that needs reversing. "Good bond" is already OK.
	BondPtr badBond = _trappedToCTerminus;
	BondPtr goodBond = _trappedToNTerminus;

	// If it's flipped we want the opposite.
	if (_flipNTerminus)
	{
		badBond = _trappedToNTerminus;
		goodBond = _trappedToCTerminus;
	}

	// Store the position of the new anchor atom (N).
	_anchorStart = goodBond->getStaticPosition();
	_majorStart = badBond->getStaticPosition();
	// Get position of heavy start before we flip the bond.

	// Perhaps we chose a stupid anchor?
	if (badBond->downstreamAtomGroupCount() != 1 || badBond->downstreamAtomCount(0) == 0)
	{
		shout_at_user("Cannot anchor a bond which has no/2+ downstream groups");
	}

	// During reversal of the bad bond, we need to store the info of the
	// next bond downstream prior to reversal.
	_newHeavyStoredGroup = ToBondPtr(badBond->downstreamAtom(0, 0)->getModel())->getBondGroup(0);

	// we can now flip the bad bond with impunity.
	badBond->reverse(BondPtr());

	// Remove all the downstream atoms, they need to be re-calculated.
	badBond->getBondGroup(0)->atoms.clear();

	/* Add all the upstream bonds as downstream bonds during semi-manual reversal.
	 * Primary atom needs to be the major atom to continue main chain. */
	badBond->addDownstreamAtom(goodBond->getMajor(), 0);
	_minor = goodBond->getMinor();

	/* Add the rest in reverse order to maintain chirality */
	for (int j = goodBond->downstreamAtomCount(0) - 1; j > 0; j--)
	{
		AtomPtr atom = goodBond->downstreamAtom(0, j);
		badBond->addDownstreamAtom(atom, 0);
	}

	/* One last sanity check for good measure */
	if (_trappedToNTerminus->getMinor() != _trappedToCTerminus->getMinor())
	{
		std::cout << "Oops, it's " << _trappedToNTerminus->shortDesc()
		<< " and " << _trappedToCTerminus->shortDesc() << std::endl;
	}
}

void Anchor::activate()
{
	// This tells the bonds not to respond to propagation of changes anymore.
	_trappedToCTerminus->setAnchored();
	_trappedToNTerminus->setAnchored();

	/* As in the prior function */
	BondPtr badBond = _trappedToCTerminus;
	BondPtr goodBond = _trappedToNTerminus;

	if (_flipNTerminus)
	{
		badBond = _trappedToNTerminus;
		goodBond = _trappedToCTerminus;
	}

	// We haven't error checked goodBond yet! */

	BondGroup *badGroup = badBond->getBondGroup(0);
	BondGroup *goodGroup = goodBond->getBondGroup(0);

	if (badBond->downstreamAtomGroupCount() == 0)
	{
		shout_at_user("Cannot anchor a bond which has no downstream groups");
	}

	/* We need to reverse the stored information in the bad bond */

	vec3 lPos = goodGroup->staticSample[0].old_start;
	vec3 miPos = goodGroup->staticSample[0].start;
	vec3 maPos = badGroup->staticSample[0].start;
	vec3 hPos = _newHeavyStoredGroup->staticSample[0].start;

	/* Bad group static position fixed. NB: absolute pos?? */
	BondSample *staticSample = &badGroup->staticSample[0];
	staticSample->start = miPos;
	staticSample->old_start = maPos;
	staticSample->basis = badBond->makeTorsionBasis(hPos, maPos, miPos, lPos,
													&staticSample->torsion);

	/* Switch over the major and minor positions, and re-calculate torsion
	 * bases from scratch for each stored sample in the bad bonds. */
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

		if (sample->torsion != sample->torsion)
		{
			std::cout << "!!" << std::endl;
		}

		sample->old_start = sample->start;
		sample->start = goodSample->start;

	}

	/* Finally, allow the anchored nitrogen to use the Anchor as its model. */
	AtomPtr nitrogen = _trappedToNTerminus->getMinor();
	nitrogen->setModel(shared_from_this());
}

BondPtr Anchor::getAppropriateBond(bool reverse)
{
	if (!_callingBond)
	{
		return _trappedToNTerminus;
	}
	else if (!reverse && _trappedToCTerminus->getMajor() == _callingBond->getMinor())
	{
		return _trappedToNTerminus;
	}
	else if (!reverse && _trappedToNTerminus->getMajor() == _callingBond->getMinor())
	{
		return _trappedToCTerminus;
	}
	else if (reverse && _trappedToCTerminus->getMajor() == _callingBond->getMajor())
	{
		return _trappedToNTerminus;
	}
	else if (reverse && _trappedToNTerminus->getMajor() == _callingBond->getMajor())
	{
		return _trappedToCTerminus;
	}
	else
	{
		return _trappedToNTerminus;
		shout_at_helen("Calling bond is not appropriate!");
		return BondPtr();
	}
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
