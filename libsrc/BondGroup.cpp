//
//  BondGroup.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "BondGroup.h"
#include "Bond.h"

void BondGroup::addProperties()
{
	for (size_t i = 0; i < bondCount(); i++)
	{
		addReference("bond", bondPtr(i));
	}
	
	addIntProperty("group", &_group);
}

void BondGroup::linkReference(BaseParserPtr object, std::string category)
{
	if (category == "bond")
	{
		BondPtr bond = ToBondPtr(object);
		addBond(bond);
	}
}

BondPtr BondGroup::bondPtr(int i)
{
	return _bonds[i]->shared_from_this();
}


