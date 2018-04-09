// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#include "Thing.h"
#include "VScript.h"
#include "Parser.h"
#include "FileReader.h"
#include <sstream>

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
		std::cout << "Class " << type << " doesn't exist." << std::endl;
		throw VErrorInvalidThingType;
	}
	
	_type = ThingParser;
}

std::string Thing::description()
{
	std::ostringstream stream;	

	switch (_type)
	{
		case ThingDouble:
		stream << _doubleValue;
		break;
		
		case ThingString:
		stream << _stringValue;
		break;
		
		case ThingInt:
		stream << _intValue;
		break;
		
		case ThingParser:
		stream << _parser->getAbsolutePath();		
		break;
		
		default:
		break;
	}
	
	return stream.str();
}

void Thing::printDescription()
{
	std::cout << description() << std::endl;
}

std::vector<std::string> splitContents(std::string contents)
{
	std::vector<std::string> params = split(contents, ',');	
	
	for (int i = 0; i < params.size(); i++)
	{
		trim(params[i]);
	}
	
	return params;
}

/* This will need to be moved into Thing */
ThingPtr Thing::dealWithFunction(std::string function, std::string contents)
{
	std::vector<std::string> params = splitContents(contents);
	
	/* Some special function calls... get_molecule, or count_molecule for
 * 	* example - we can try and find it now. */
	
	bool get = false;
	bool count = false;
	char *pos = &function[0];
	
	if (function.compare(0, 4, "get_") == 0)
	{
		get = true;
		pos += 4;
	}
	else if (function.compare(0, 6, "count_") == 0)
	{
		count = true;
		pos += 6;
	}
	
	if (_type != ThingParser && (get || count))
	{
		throw VErrorSimpleTypeFunctionCall;
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
			if (params.size() < 1)
			{
				throw VErrorMissingParameter;	
			}
			
			int get_num = atoi(params[0].c_str());
			
			if (get_num >= childCount)
			{
				throw VErrorBeyondArrayBounds;
			}
			
			ParserPtr parser = getParserValue()->getChild(object, get_num);
			right->setThingType(ThingParser);
			right->setParserValue(parser);
			return right;
		}

		if (dPtr) right->setDoubleValue(*dPtr);
		if (iPtr) right->setIntValue(*iPtr);
		if (sPtr) right->setStringValue(*sPtr);
		
		return right;
	}
	
	/* It was neither a getter nor a counter */
	
	if (function == "print")
	{
		printDescription();
		return ThingPtr();
	}
	
	std::cout << "Could not identify function?" << std::endl;
	throw VErrorMissingImplementation;
	
	
	return ThingPtr();
}

void Thing::addThing(ThingPtr right)
{
	if (this->getThingType() == ThingParser)
	{
		throw VErrorImmutableObject;
	}

	if (this->getThingType() == ThingString)
	{
		std::string addition = right->description();
		_stringValue += addition;
		return;
	}
	
	if (this->isNumber() && !right->isNumber())
	{
		throw VErrorInappropriateOperation;	
	}
	
	if (this->isNumber() && right->isNumber())
	{
		if (getThingType() == ThingInt)
		{
			_intValue += right->getIntValue();
		}

		if (getThingType() == ThingDouble)
		{
			_doubleValue += right->getDoubleValue();
		}
		
		return;
	}
	
	throw VErrorInappropriateOperation;
}



