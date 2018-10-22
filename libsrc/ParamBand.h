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

/**
 * \class ParamBand
 * \brief An object which looks after a related set of parameters for
 * refinement.
 *
 * A band of parameters all of the same type can be stored and changed as one
 * using this parameter band.
 * Initialise, and set the private getter/setter for the parameter type.
 * Add each object for which the getter/setter pertains.
 * Do not forget to "prepare" the object, which will take local copies
 * of the original values.
 * You can also use this to store parameters of a given type and then
 * reapply them by using setGlobalParam(this, 0), instead of using them
 * for refinement.
 */

class ParamBand
{
public:
	ParamBand();
	ParamBand(ParamBand &other);
	
	void setPrivateGetter(Getter getter)
	{
		_getter = getter;
	}
	
	void setPrivateSetter(Setter setter)
	{
		_setter = setter;
	}
	
	void addObject(void *object, double weight);
	
	/** Adds 'value' to each parameter in the object list. */
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
	
	void *object(int i)
	{
		return _objects[i].object;
	}
	
	double baseValueForObject(int i)
	{
		return _objects[i].baseValue;
	}
	
	void prepare();
private:
	Getter _getter;
	Setter _setter;
	double _paramValue;

	std::vector<ObjectWeight> _objects;
};

#endif
