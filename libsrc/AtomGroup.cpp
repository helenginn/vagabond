//
//  AtomGroup.cpp
//  vagabond
//
//  Created by Helen Ginn on 25/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "AtomGroup.h"
#include "Atom.h"
#include "Element.h"
#include "Bond.h"
#include <sstream>

AtomPtr AtomGroup::findAtom(std::string atomType)
{
	for (int i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getAtomName() == atomType)
		{
			return atom(i);
		}
	}

	return AtomPtr();
}

AtomList AtomGroup::findAtoms(std::string atomType)
{
	AtomList list;

	for (int i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getAtomName() == atomType)
		{
			list.push_back(atom(i));
		}
	}

	return list;
}

double AtomGroup::totalElectrons()
{
	double total = 0;

	for (int i = 0; i < atomCount(); i++)
	{
		total += atom(i)->getElement()->electronCount();
	}

	return total;
}

std::string AtomGroup::getPDBContribution()
{
	std::ostringstream stream;

	for (int i = 0; i < bondCount(); i++)
	{
		BondPtr aBond = bond(i);
		stream << aBond->getPDBContribution();
	}

	return stream.str();
}

void AtomGroup::setUseAbsolute()
{
	for (int i = 0; i < atomCount(); i++)
	{
		atom(i)->setKeepModel();
	}
}


AtomGroup::AtomGroup()
{
	_beenTied = false;
	_timesRefined = 0;
}