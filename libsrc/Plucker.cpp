//
//  Plucker.cpp
//  vagabond
//
//  Created by Helen Ginn on 05/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "Plucker.h"
#include <math.h>
#include <iostream>

Plucker::Plucker()
{
	_granularity = 0.1;
}


void Plucker::_addPluckable(void *pluck, double weight)
{
	int chances = lrint(weight / _granularity);

	for (int i = 0; i < chances; i++)
	{
		_pluckable.push_back(pluck);
	}
}

void Plucker::printSummary()
{
	std::cout << "Plucker has " << _pluckable.size()
	<< " choices to choose from." << std::endl;
}




