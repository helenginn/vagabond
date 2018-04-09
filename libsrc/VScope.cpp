// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#include "VScope.h"
#include "VScript.h"
#include "LeftThing.h"
#include <iostream>

VScope::VScope()
{

}


void VScope::addLeftThing(LeftThingPtr thing)
{
	std::string name = thing->getThingName();
	
	if (_things.count(name) > 0)
	{
		throw VErrorThingRedeclaration;
	}
	
	_things[name] = thing;
}

LeftThingPtr VScope::findThing(std::string name)
{
	if (_things.count(name) == 0)
	{
		return LeftThingPtr();
	}
	
	return _things[name];
}
