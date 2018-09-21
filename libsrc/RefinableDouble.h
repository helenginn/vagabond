//
//    RefinableDouble.h
//    Vagabond - stuff

//    Copyright (C) 2018  Helen Ginn
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __vagabond__RefinableDouble__
#define __vagabond__RefinableDouble__

#include <vector>
#include <iostream>
#include "RefinementStrategy.h"

typedef double (*TaggedGetter)(void *, int);

class RefinableDouble
{
public:
	RefinableDouble(int num, Clusterable *parent)
	{
		_value = 0;
		_num = num;
		_parent = parent;
		_gradient = NULL;
	}
	
	void setGradientFunction(TaggedGetter gradient)
	{
		_gradient = gradient;
	}

	static double getDouble(void *object)
	{
		RefinableDouble *val = static_cast<RefinableDouble *>(object);
		return val->_value;
	}

	static double getGradient(void *object)
	{
		RefinableDouble *val = static_cast<RefinableDouble *>(object);
		return (*(val->_gradient))(val->_parent, val->_num);
	}

	static void setDouble(void *object, double value)
	{
		RefinableDouble *val = static_cast<RefinableDouble *>(object);
		val->_value = value;
	}

private:
	Clusterable *_parent;
	int _num;
	double _value;
	TaggedGetter _gradient;
};


#endif
