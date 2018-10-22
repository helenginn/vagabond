//
//  BondGroup.h
//  vagabond
//
//  Created by Helen Ginn on 23/07/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

/**
* \class BondGroup
*
* \brief The BondGroup class holds an array of downstream bonds for a given 
*  bond. One parent bond will have one BondGroup for every alternative conformer.
*
*/

#ifndef __vagabond__BondGroup__
#define __vagabond__BondGroup__

#include <stdio.h>
#include <string>
#include <vector>
#include "Parser.h"
#include "FileReader.h"

class BondGroup : public Parser
{
public:
	BondGroup(int group)
	{
		_group = group;
	}

	BondGroup()
	{
		_group = -1;
	}
	
	void addBond(BondPtr bond)
	{
		_bonds.push_back(&*bond);
	}
	
	void addBond(Bond *bond)
	{
		_bonds.push_back(bond);
	}
	
	size_t bondCount()
	{
		return _bonds.size();
	}
	
	BondPtr bondPtr(int i);

	Bond *bond(int i)
	{
		return _bonds[i];
	}
	
	void removeBond(int i)
	{
		_bonds.erase(_bonds.begin() + i);
	}

	void clearBonds()
	{
		_bonds.clear();
	}
protected:
	virtual void addProperties();
	virtual void linkReference(ParserPtr object, std::string category);
	virtual std::string getParserIdentifier()
	{
		return "group_" + i_to_str(_group);
	}
	
	virtual std::string getClassName()
	{
		return "BondGroup";
	}

private:
	int _group;	
	std::vector<Bond *> _bonds;
};

#endif
