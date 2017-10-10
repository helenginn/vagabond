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
#include "Absolute.h"
#include "shared_ptrs.h"
#include <vector>

class Anchor : public Absolute
{
public:
	Anchor(BondPtr inheritDownstreamBond, BondPtr inheritParentBond);

	void setCallingBond(Bond *call)
	{
		_callingBond = call;
	}

	virtual vec3 getStaticPosition()
	{
		return _staticPosition;
	}

	virtual vec3 getAbsolutePosition()
	{
		return _absPosition;
	}

	virtual std::vector<BondSample> *getManyPositions(BondSampleStyle style)
	{
		return &_manyPositions;
	}

	virtual std::string getClassName()
	{
		return "Anchor";
	}
private:

	vec3 _staticPosition;
	vec3 _absPosition;
	std::vector<BondSample> _manyPositions;

	/* Anchor will copy each of the input bonds and take a static
	 * copy of them. Anchor will then behave like a bond, but return
	 * the appropriate copy depending on what _callingBond is */
	Bond *_trappedToNTerminus;
	Bond *_trappedToCTerminus;

	Bond *_callingBond;
};

#endif /* defined(__vagabond__Anchor__) */
