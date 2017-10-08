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

	for (int i = 0; i < atomCount(); i++)
	{
		stream << atom(i)->getPDBContribution();
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

void AtomGroup::addAtomsFrom(AtomGroupPtr child)
{
	for (int i = 0; i < child->atomCount(); i++)
	{
		addAtom(child->atom(i));
	}
}

double AtomGroup::getAverageBFactor(bool initial)
{
	double sum = 0;
	double count = 0;

	for (int i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getElement()->electronCount() <= 1)
		{
			continue;
		}

		if (initial)
		{
			sum += atom(i)->getInitialBFactor();
			count++;
		}
		else
		{
			if (atom(i)->getModel()->isBond())
			{
				BondPtr bond = ToBondPtr(atom(i)->getModel());
				double val = bond->getMeanSquareDeviation();
				sum += val;
				count++;
			}
		}
	}

	return sum / count;
}

AtomGroup::AtomGroup()
{
	_beenTied = false;
	_timesRefined = 0;
}
