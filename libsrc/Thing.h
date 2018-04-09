// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#ifndef __vagabond__Thing__
#define __vagabond__Thing__

#include "shared_ptrs.h"

typedef enum
{
	ThingUnassigned,
	ThingParser,
	ThingInt,
	ThingString,
	ThingDouble,
	ThingInvalid,
} ThingType;

class Thing
{
public:
	Thing();
	
	void setThingType(std::string type);
	void setThingType(ThingType type)
	{
		_type = type;	
	}
	
	ThingType getThingType()
	{
		return _type;	
	}
	
	void setDoubleValue(double value)
	{
		_doubleValue = value;
	}

	void setParserValue(ParserPtr value)
	{
		_parser = value;
	}

	
	void setIntValue(int value)
	{
		_intValue = value;
	}

	void setStringValue(std::string value)
	{
		_stringValue = value;
	}
	
	int getIntValue()
	{
		return _intValue;
	}
	
	double getDoubleValue()
	{
		return _doubleValue;
	}
	
	std::string getStringValue()
	{
		return _stringValue;
	}
	
	ParserPtr getParserValue()
	{
		return _parser;
	}
protected:
	ThingType _type;
	
private:	
	ParserPtr _parser;

	int _intValue;
	double _doubleValue;
	std::string _stringValue;
};

#endif
