//
//  FlexGlobal.hpp
//  vagabond
//
//  Created by Helen Ginn on 27/12/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#ifndef FlexGlobal_hpp
#define FlexGlobal_hpp

#include <stdio.h>
#include "shared_ptrs.h"

class FlexGlobal
{
public:
	FlexGlobal();

	void setAtomGroup(AtomGroupPtr group)
	{
		_atomGroup = group;
	}

	static double score(void *object);

	void maximiseIsotropy();
	void setTargetBFactor(double value)
	{
		_targetIsoB = value;
	}
private:
	double _targetIsoB;

	double notStaticScore();
	AtomGroupPtr _atomGroup;
};

#endif /* FlexGlobal_hpp */
