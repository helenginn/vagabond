// Vagabond : bond-based macromolecular model refinement
// Copyright (C) 2017-2018 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#include "ParamBand.h"
#include "Shouter.h"
#include <iostream>

ParamBand::ParamBand()
{
	_paramValue = 0;
	_allWeight = 1;
	_equalDiv = false;
	_min = -FLT_MAX;
	_max = +FLT_MAX;
}

ParamBand::ParamBand(ParamBand &other)
{
	_getter = other._getter;
	_setter = other._setter;
	_paramValue = other._paramValue;
	_objects = other._objects;
	_equalDiv = false;
	_allWeight = other._allWeight;
	_min = other._min;
	_max = other._max;
}

void ParamBand::addObject(void *object, double weight)
{
	ObjectWeight pair;
	pair.object = object;
	pair.weight = weight;
	pair.baseValue = 0;
	_objects.push_back(pair);
}

void ParamBand::prepare()
{
	if (_getter == NULL)
	{
		shout_at_helen("Need to set global getter for param band");
	}
	
	for (int i = 0; i < objectCount(); i++)
	{
		ObjectWeight obj = _objects[i];
		double value = (*_getter)(obj.object);
		_objects[i].baseValue = value;
	}
}

void ParamBand::setGlobalParam(void *object, double value)
{
	ParamBand *me = static_cast<ParamBand *>(object);
	double divValue = value / (double)me->objectCount();
	
	if (me->_equalDiv)
	{
		divValue = value;
	}
	
	for (int i = 0; i < me->objectCount(); i++)
	{
		double base = me->_objects[i].baseValue;
		double newVal = (base + divValue) * me->_allWeight;
		
		if (newVal > me->_max)
		{
			newVal = me->_max;
		}
		
		if (newVal < me->_min)
		{
			newVal = me->_min;
		}

		(*me->_setter)(me->_objects[i].object, newVal);
	}
	
	me->_paramValue = value;
}

double ParamBand::getGlobalParam(void *object)
{
	ParamBand *me = static_cast<ParamBand *>(object);
	return me->_paramValue;
}
