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

#ifndef __vagabond__ParamBand_h__
#define __vagabond__ParamBand_h__

#include "RefinementStrategy.h"

typedef struct
{
	void *object;
	double weight;
	double baseValue;
} ObjectWeight;

class ParamBand
{
public:
	ParamBand();
	
	void setPrivateGetter(Getter getter)
	{
		_getter = getter;
	}
	
	void setPrivateSetter(Setter setter)
	{
		_setter = setter;
	}
	
	void addObject(void *object, double weight);
	
	static void setGlobalParam(void *object, double value);
	static double getGlobalParam(void *object);
	
	static bool has_more_objects(ParamBandPtr &one, ParamBandPtr &two)
	{
		return (one->objectCount() > two->objectCount());
	}
	
	size_t objectCount()
	{
		return _objects.size();
	}
	
	void prepare();
private:
	Getter _getter;
	Setter _setter;
	double _paramValue;

	std::vector<ObjectWeight> _objects;
};

#endif
