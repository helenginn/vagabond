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

AtomPtr AtomGroup::findAtom(std::string atomType, std::string confID)
{
	AtomList atoms = findAtoms(atomType);

	for (int i = 0; i < atoms.size(); i++)
	{
		if (atoms[i].expired())
		{
			continue;
		}

		if (atoms[i].lock()->getAlternativeConformer() == confID)
		{
			return atoms[i].lock();
		}
	}

	return AtomPtr();
}

std::map<std::string, int> AtomGroup::conformerMap()
{
	std::map<std::string, int> conformerList;

	for (int i = 0; i < atomCount(); i++)
	{
		std::string conformer = atom(i)->getAlternativeConformer();

		if (!conformerList.count(conformer))
		{
			conformerList[conformer] = 0;
		}

		conformerList[conformer]++;
	}

	return conformerList;
}

int AtomGroup::conformerCount()
{
	std::map<std::string, int> conformerList = conformerMap();

	return conformerList.size();
}

std::string AtomGroup::conformer(int i)
{
	std::map<std::string, int> conformerList = conformerMap();
	std::map<std::string, int>::iterator it = conformerList.begin();

	for (int j = 0; j < i; j++) it++;

	return it->first;
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

std::string AtomGroup::getPDBContribution(PDBType pdbType, CrystalPtr crystal)
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

				if (atom(i)->getWeighting() <= 0)
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
		if (crystal)
		{
			stream << atom(i)->anisouPDBLine(crystal);
		}
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
			ToBondPtr(atom(i)->getModel())->calculateInitialMagicAxis();
		}
	}
}

AtomList AtomGroup::topLevelAtoms()
{
	if (!atomCount()) return AtomList();

	AtomPtr topAtom = atom(0);

	while (true)
	{
		if (!topAtom->getModel()->isBond())
		{
			break;
		}

		BondPtr bond = ToBondPtr(topAtom->getModel());

		if (!hasAtom(bond->getMajor()))
		{
			break;
		}

		topAtom = bond->getMajor();
	}

	AtomList list;
	list.push_back(topAtom);

	return list;
}

bool AtomGroup::hasAtom(AtomPtr anAtom)
{
	bool found = false;

	if (!anAtom) return false;

	for (int i = 0; i < atomCount(); i++)
	{
		if (atom(i) == anAtom)
		{
			found = true;
		}
	}

	return found;
}

void AtomGroup::refine(CrystalPtr target, RefinementType rType)
{
	AtomList topAtoms = topLevelAtoms();

	ScoreType scoreType = ScoreTypeModelPos;
	int maxTries = 60;
	int bondNum = 4;

	if (rType == RefinementFine)
	{
		scoreType = ScoreTypeMultiply;
		maxTries = 1;
		bondNum = 5;
	}

	for (int n = 0; n < topAtoms.size(); n++)
	{
		AtomPtr topAtom = topAtoms[n].lock();

		if (n == 1) std::cout << "'" << std::flush;

		while (hasAtom(topAtom))
		{
			if (!topAtom->getModel()->isBond())
			{
				break;
			}

			BondPtr bond = ToBondPtr(topAtom->getModel());

			int groups = bond->downstreamAtomGroupCount();

			if (!groups)
			{
				break;
			}

			if (!bond->isRefinable())
			{
				break;
			}

			for (int k = 0; k < 1; k++)
			{
				if (shouldRefineMagicAxis(bond))
				{
					bond->calculateMagicAxis();
				}

				bool changed = true;
				int count = 0;

				BondPtr topBond;
				while (changed && count < maxTries)
				{
					bond->setActiveGroup(k);
					setupNelderMead();
					setCrystal(target);
					topBond = setupTorsionSet(bond, k, bondNum,
											  deg2rad(4), deg2rad(0.04));
					setScoreType(scoreType);

					setSilent();

					setJobName("torsion_" +  bond->shortDesc());
					changed = sample();
					count++;
				}

				if (!topBond)
				{
					topAtom = AtomPtr();
					continue;
				}

				topAtom = topBond->getMinor();
			}
		}
	}
}
