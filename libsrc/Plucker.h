//
//  Plucker.h
//  vagabond
//
//  Created by Helen Ginn on 05/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#ifndef __vagabond__Plucker__
#define __vagabond__Plucker__

#include <vector>
#include <iostream>
#include <math.h>

class Plucker
{
public:
	Plucker();
	
	/** Each weight for Plucked objects will be divided by this number, which
	* defines the number of chances the object has of being plucked relative
	* to the others. */
	void setGranularity(double gran)
	{
		_granularity = gran;
	}

	/** Add an object with an appropriate weighting to the list of pluckable
	* objects */
	template <class Plucked> void addPluckable(Plucked *pluck, double weight)
	{
		_addPluckable((void *)pluck, weight);
	}
	
	size_t pluckCount()
	{
		return _pluckable.size();
	}
	
	/** Grab an object from the weighted chances and cast to desired return
	* type */
	void *pluck()
	{
		if (_pluckable.size() == 0)
		{
			return NULL;
		}
		
		int random = rand() % (_pluckable.size());

		return _pluckable[random];
	}
	
	void printSummary();
private:	
	void _addPluckable(void *pluck, double weight);
	std::vector<void *> _pluckable;

	double _granularity;
};

#endif
