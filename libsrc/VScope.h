// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#include <map>

#include "shared_ptrs.h"

/** Looks after a bunch of LeftThings */

class VScope
{
public:
	VScope();
	
	void addLeftThing(LeftThingPtr thing);
	LeftThingPtr findThing(std::string name);
private:	
	
	std::map<std::string, LeftThingPtr> _things;
};
