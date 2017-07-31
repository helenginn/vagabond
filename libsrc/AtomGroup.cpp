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

double AtomGroup::totalElectrons()
{
	double total = 0;

	for (int i = 0; i < atomCount(); i++)
	{
		total += atom(i)->getElement()->electronCount();
	}

	return total;
}

void AtomGroup::getPDBContribution()
{
	for (int i = 0; i < atomCount(); i++)
	{
		ModelPtr model = atom(i)->getModel();

		if (model->getClassName() == "Bond")
		{
			BondPtr bond = std::static_pointer_cast<Bond>(model);


			bond->getPDBContribution();
		}
	}

}
