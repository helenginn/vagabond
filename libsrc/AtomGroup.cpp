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
#include <iomanip>

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

std::string AtomGroup::getPDBContribution(PDBType pdbType)
{
	std::ostringstream stream;
	int numConf = 0;

	if (!atomCount())
	{
		return "";
	}

	if (pdbType == PDBTypeEnsemble)
	{
		/* Get the total number of conformers to worry about */
		std::vector<BondSample> *samples = atom(0)->getModel()->getManyPositions(BondSampleThorough);

		numConf = samples->size();

		for (int j = 0; j < numConf; j++)
		{
			stream << "MODEL " << std::setw(8) << j << std::setw(66) << " " << std::endl;

			for (int i = 0; i < atomCount(); i++)
			{
				if (!atom(i)->getMonomer())
				{
					continue;
				}

				stream << atom(i)->getPDBContribution(j);
			}
			
			stream << "TER" << std::setw(80) << " " << std::endl;
			stream << "ENDMDL" << std::setw(80) << " " << std::endl;
		}

		return stream.str();
	}

	for (int i = 0; i < atomCount(); i++)
	{
		bool samePos = (pdbType == PDBTypeSamePosition);
		bool sameB = (pdbType == PDBTypeSameBFactor);
		stream << atom(i)->averagePDBContribution(samePos, sameB);
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

double AtomGroup::getAverageDisplacement()
{
	double sum = 0;
	double count = 0;

	for (int i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getElement()->electronCount() <= 1)
		{
			continue;
		}

		double val = atom(i)->posDisplacement();

		sum += val;
		count++;
	}

	return sum / count;
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
}

void AtomGroup::propagateChange()
{
	for (int i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getModel()->isBond())
		{
			ToBondPtr(atom(i)->getModel())->propagateChange();
		}
	}
}

int AtomGroup::totalElectrons(int *fcWeighted)
{
	double sum = 0;
	double weighted = 0;

	for (int i = 0; i < atomCount(); i++)
	{
		double e = atom(i)->getElement()->electronCount();
		sum += e;
		double weight = atom(i)->getWeighting();
		weighted += e * weight;
	}

	*fcWeighted = weighted;

	return sum;
}

void AtomGroup::setWeighting(double value)
{
	for (int i = 0; i < atomCount(); i++)
	{
		atom(i)->setWeighting(value);
	}
}

void AtomGroup::resetMagicAxes()
{
	for (int i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getModel()->isBond())
		{
			ToBondPtr(atom(i)->getModel())->calculateMagicAxis();
		}
	}
}
