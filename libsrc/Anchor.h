//
//  Anchor.h
//  vagabond
//
//  Created by Helen Ginn on 11/09/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Anchor__
#define __vagabond__Anchor__

#include <stdio.h>
#include "Bond.h"
#include "shared_ptrs.h"
#include <vector>

class Anchor : public Bond
{
public:
	Anchor(BondPtr inheritDownstreamBond, BondPtr inheritParentBond);

	void setCallingBond(Bond *call)
	{
		_callingBond = call;
	}

	virtual vec3 getStaticPosition()
	{
		return getAppropriateBond()->getStaticPosition();
	}

	virtual vec3 getAbsolutePosition()
	{
		return getAppropriateBond()->getAbsolutePosition();
	}

	virtual std::vector<BondSample> *getManyPositions(BondSampleStyle style)
	{
		return getAppropriateBond()->getManyPositions(style);
	}

	virtual std::string getClassName()
	{
		return "Anchor";
	}

	/* Which trapped bond is required for the _callingBond? */
	/* If the calling bond happens to be reversed, pass reverse = true */
	BondPtr getAppropriateBond(bool reverse = false);
	void activate();

	static BondPtr sanitiseBond(Bond *myself, BondPtr model);
private:

	/* Anchor will copy each of the input bonds and take a static
	 * copy of them. Anchor will then behave like a bond, but return
	 * the appropriate copy depending on what _callingBond is */
	BondPtr _trappedToNTerminus;
	BondPtr _trappedToCTerminus;
	bool _flipNTerminus;

	/* Temporary holding of static start position for anchor, heavy atom and
	 * bond group */
	vec3 _anchorStart;
	vec3 _majorStart;
	BondGroup *_newHeavyStoredGroup;

	Bond *_callingBond;
};

#endif /* defined(__vagabond__Anchor__) */
