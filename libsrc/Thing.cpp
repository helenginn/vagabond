// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#include "Thing.h"
#include "VScript.h"
#include "Parser.h"

Thing::Thing()
{
	_type = ThingUnassigned;
}

void Thing::setThingType(std::string type)
{
	/* A simple type */
	if (type == "int")
	{
		_type = ThingInt;
		return;
	}
	else if (type == "double")
	{
		_type = ThingDouble;
		return;	
	}
	else if (type == "string")
	{
		_type = ThingString;
		return;
	}

	int count = Parser::classCount(type);
	
	if (!count)
	{
		throw VErrorInvalidThingType;
	}
	
	_type = ThingParser;
}

