//
//  FlexRegion.hpp
//  vagabond
//
//  Created by Helen Ginn on 04/01/2018.
//  Copyright © 2018 Strubi. All rights reserved.
//

#ifndef FlexRegion_hpp
#define FlexRegion_hpp

#include <stdio.h>
#include "shared_ptrs.h"

class FlexRegion
{
public:
	FlexRegion();

	void addBond(BondPtr bond, int prevBondCount);
	void addSingleBondParameters();
	void setup();
	void sample();

	int bondCount()
	{
		return _bonds.size();
	}
private:
	void addSingleBondParameter(int i);
	static double score(void *object)
	{
		return static_cast<FlexRegion *>(object)->getScore();
	}

	double getScore();

	RefinementStrategyPtr _strategy;

	std::vector<BondPtr> _bonds;
};

#endif /* FlexRegion_hpp */
