// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#ifndef __vagabond__Thing__
#define __vagabond__Thing__

#include "shared_ptrs.h"
#include "VScript.h"

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
	
	ThingPtr dealWithFunction(std::string function, 
	                          std::vector<ThingPtr> things);
	void printDescription();
	std::string description();
	
	virtual void addThing(ThingPtr right);
	
	VScriptComparison compareToThing(ThingPtr other);
	void setThingType(std::string type);
	void setThingType(ThingType type)
	{
		_type = type;	
	}
	
	bool isNumber()
	{
		return _type == ThingInt || _type == ThingDouble;
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
		if (_type == ThingDouble)
		{
			return _doubleValue;
		}
		
		return _intValue;
	}
	
	double getDoubleValue()
	{
		if (_type == ThingInt)
		{
			return _intValue;
		}
		
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
	int _intValue;
	double _doubleValue;
	std::string _stringValue;
	
private:	
	ParserPtr _parser;

};

template <class T>
VScriptComparison compare(T &one, T &two)
{
	if (one == two)
	{
		return VCompEqual;
	}
	else if (one > two)
	{
		return VCompGreaterThan;
	}
	else if (one < two)
	{
		return VCompLessThan;
	}
}


#endif
