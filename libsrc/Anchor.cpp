//
//  Anchor.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/09/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Anchor.h"
#include "Bond.h"

Anchor::Anchor(BondPtr inheritDownstreamBond, BondPtr inheritParentBond)
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
	 * The position of N is described by downstream position.
	 * the position of CA is described by parent position.
	 */

	inheritDownstreamBond->getDistribution(); // update mean position
	inheritParentBond->getDistribution(); // update mean position

	/* Get the distributions of all positions from both parent and downstream bond */
	std::vector<BondSample> *manyParentPositions = inheritParentBond->getManyPositions(BondSampleThorough);
	std::vector<BondSample> *manyDownstreamPositions = inheritDownstreamBond->getManyPositions(BondSampleThorough);

	/* Set our internal static position to the parent bond static position. */
	_staticPosition = inheritDownstreamBond->getStaticPosition();

	/* Set our absolute static position to the parent bond static position.
	 * Will be returned when asking for PDB, atom position etc. */
	_absPosition = inheritDownstreamBond->getStaticPosition();

	ModelPtr grandParent = inheritParentBond->downstreamAtom(0, 0)->getModel();
	BondPtr grandBond = ToBondPtr(grandParent);
	grandBond->getDistribution();
	std::vector<BondSample> *manyGrandPositions = grandBond->getManyPositions(BondSampleThorough);

	for (int i = 0; i < manyDownstreamPositions->size(); i++)
	{
		BondSample parentSample = manyParentPositions->at(i);
		parentSample.basis.vals[2] = -parentSample.basis.vals[2];
		parentSample.basis.vals[5] = -parentSample.basis.vals[5];
		parentSample.basis.vals[8] = -parentSample.basis.vals[8];
		parentSample.basis.vals[1] = -parentSample.basis.vals[1];
		parentSample.basis.vals[4] = -parentSample.basis.vals[4];
		parentSample.basis.vals[7] = -parentSample.basis.vals[7];
		parentSample.start = parentSample.old_start;
		parentSample.old_start = manyGrandPositions->at(i).start;
		_manyPositions.push_back(parentSample);
	}
}