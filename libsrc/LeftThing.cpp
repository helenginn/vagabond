// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#include "LeftThing.h"
#include "Thing.h"
#include <string>
#include "Parser.h"
#include "VScript.h"

void LeftThing::setThingName(std::string name)
{
	_name = name;
}


void LeftThing::setThing(ThingPtr right)
{
	checkTypesAgree(right);
	
	if (_type == ThingInt)
	{
		int val = right->getIntValue();
		setIntValue(val);
	}
	else if (_type == ThingDouble)
	{
		double val = right->getDoubleValue();
		setDoubleValue(val);
	}
	else if (_type == ThingString)
	{
		std::string val = right->getStringValue();
		setStringValue(val);
	}
	else if (_type == ThingParser)
	{
		ParserPtr val = right->getParserValue();
		setParserValue(val);
	}
}

void LeftThing::checkTypesAgree(ThingPtr right)
{
	if (isNumber() && right->isNumber())
	{
		return;
	}
	
	if (_type != right->getThingType())
	{
		throw VErrorTypeMismatch;
	}
}

