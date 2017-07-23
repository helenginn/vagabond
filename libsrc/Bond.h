//
//  Bond.h
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Bond__
#define __vagabond__Bond__

#include <stdio.h>
#include "Model.h"
#include <vector>

class Bond : public Model
{
public:
	Bond(AtomPtr major, AtomPtr minor);
	void setup();

	AtomPtr getMajor()
	{
		return _major.lock();
	}

	AtomPtr getMinor()
	{
		return _minor.lock();
	}
	
	virtual FFTPtr getDistribution();
public:
	static double getVoxelValue(void *obj, double x, double y, double z);

private:

	AtomWkr _major;
	AtomWkr _minor;

//	std::vector<AtomWkr> _dependencies;
};

#endif /* defined(__vagabond__Bond__) */
