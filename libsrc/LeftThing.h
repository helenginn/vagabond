// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#ifndef __vagabond__LeftThing__
#define __vagabond__LeftThing__

#include "Thing.h"

/** Left Things (in VScript), i.e. variables that the user can declare. */
class LeftThing : public Thing
{
public:
	void setThingName(std::string name);

	std::string getThingName()
	{
		return _name;
	}
	
	void setThing(ThingPtr right);
private:	
	std::string _name;

	void checkTypesAgree(ThingPtr right);
};

#endif
