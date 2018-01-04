//
//  FlexTarget.hpp
//  vagabond
//
//  Created by Helen Ginn on 27/12/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#ifndef FlexTarget_hpp
#define FlexTarget_hpp

#include <stdio.h>
#include "shared_ptrs.h"

class FlexTarget
{
public:
	FlexTarget();

	void setAtomGroup(AtomGroupPtr group)
	{
		_atomGroup = group;
	}

	static double score(void *object);

	void setTargetBFactor(double value)
	{
		_targetIsoB = value;
	}
private:
	double _targetIsoB;

	double notStaticScore();
	AtomGroupPtr _atomGroup;
};

#endif /* FlexTarget_hpp */
