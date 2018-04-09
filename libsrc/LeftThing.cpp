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

ThingPtr LeftThing::dealWithFunction(std::string function, std::string contents)
{
	if (_type != ThingParser)
	{
		throw VErrorSimpleTypeFunctionCall;
	}
	
	/* Some special function calls... get_molecule, or count_molecule for
	* example - we can try and find it now. */
	
	bool get = false;
	bool count = false;
	char *pos = &function[0];
	

	if (function.compare(0, 4, "get_"))
	{
		get = true;
		pos += 4;
	}
	else if (function.compare(0, 6, "count_"))
	{
		count = true;
		pos += 6;
	}
	
	ThingType _subThing;

	if (get || count)
	{
		/* Get the type of the thing */	
		std::string object = std::string(pos);

		/* Do we have a string, double, int or parser of this name? */
		
		double *dPtr = getParserValue()->getDoubleProperty(object);
		int *iPtr = getParserValue()->getIntProperty(object);
		std::string *sPtr = getParserValue()->getStringProperty(object);
		int childCount = getParserValue()->getChildCount(object);

		if (dPtr)
		{
			_subThing = ThingDouble;	
		}
		else if (iPtr)
		{
			_subThing = ThingInt;		
		}
		else if (sPtr)
		{
			_subThing = ThingString;
		}
		else if (childCount >= 0)
		{
			_subThing = ThingParser;
		}
		/* Well, none of those were ok */
		else if (get)
		{
			throw VErrorGetterDoesNotExist;	
		}
		else if (count)
		{
			throw VErrorCounterDoesNotExist;	
		}

		/* If we tried to count a single object, throw a separate error. */

		if (_subThing != ThingParser && count)
		{
			throw VErrorCounterInappropriate;
		}
		
		ThingPtr right = ThingPtr(new Thing());
		right->setThingType(_subThing);
		
		if (_subThing == ThingParser && count)
		{
			right->setThingType(ThingInt);
			right->setIntValue(childCount);
			return right;
		}
		else if (_subThing == ThingParser && get)
		{
			throw VErrorMissingImplementation;	
		}

		if (dPtr) right->setDoubleValue(*dPtr);
		if (iPtr) right->setIntValue(*iPtr);
		if (sPtr) right->setStringValue(*sPtr);
		
		return right;
	}
	
	/* It was neither a getter nor a counter */
	
	throw VErrorMissingImplementation;
	
	
	return ThingPtr();
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
	if (_type != right->getThingType())
	{
		throw VErrorTypeMismatch;
	}
}

