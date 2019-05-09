// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#ifndef __vagabond__StateValue__
#define __vagabond__StateValue__

#include "ParserTypes.h"
#include "shared_ptrs.h"

class StateValue
{
public:
	void addStringValue(std::string ptrName, std::string value)
	{
		_ptrName = ptrName;
		_string = value;
	}
	
	void addDoubleValue(std::string ptrName, double value)
	{
		_ptrName = ptrName;
		_double = value;
	}	
	
	void addMat3x3Value(std::string ptrName, mat3x3 value)
	{
		_ptrName = ptrName;
		_mat3x3 = value;
	}
	
	void addVec3Value(std::string ptrName, vec3 value)
	{
		_ptrName = ptrName;
		_vec3 = value;
	}

	void addVec3ArrayValue(std::string ptrName, std::vector<vec3> value)
	{
		_ptrName = ptrName;
		_vec3Array = value;
	}

	void addMat3x3ArrayValue(std::string ptrName, std::vector<mat3x3> value)
	{
		_ptrName = ptrName;
		_mat3x3Array = value;
	}
	
	void addBoolValue(std::string ptrName, bool value)
	{
		_ptrName = ptrName;
		_bool = value;
	}

	void addIntValue(std::string ptrName, int value)
	{
		_ptrName = ptrName;
		_int = value;
	}

	void applyToParser(BaseParser *parser);

private:
	void repairCustom();
	
	std::string _ptrName;

	std::string _string;
	std::string _custom;
	double _double;
	mat3x3 _mat3x3;
	vec3 _vec3;
	std::vector<vec3> _vec3Array;
	std::vector<mat3x3> _mat3x3Array;
	bool _bool;
	int _int;
};






#endif
