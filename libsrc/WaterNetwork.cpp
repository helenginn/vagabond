//
//  WaterNetwork.cpp
//  vagabond
//
//  Created by Helen Ginn on 01/07/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "WaterNetwork.h"
#include "Crystal.h"
#include "Sponge.h"
#include "Atom.h"
#include "RefinementNelderMead.h"

#define MAX_HBOND_DISTANCE 3.4

WaterNetwork::WaterNetwork()
{

}

void WaterNetwork::summary()
{
	Molecule::summary();
	std::cout << "| I am a water network with " << atomCount() 
	<< " waters." << std::endl;
}

void WaterNetwork::addProperties()
{
	Molecule::addProperties();
}

void WaterNetwork::postParseTidy()
{
	Molecule::postParseTidy();
}

bool hasWaters(AtomGroupPtr list)
{
	for (int i = 0; i < list->atomCount(); i++)
	{
		if (list->atom(i)->isHeteroAtom() &&
		    list->atom(i)->getAtomName() == "O")
		{
			return true;
		}
	}
	
	return false;
}

bool spongeToBeAdded(SpongePtr seed, std::vector<AtomPtr> &atoms,
                     std::vector<SpongePtr> &sponges)
{
	if (std::find(sponges.begin(), sponges.end(), seed) != sponges.end())
	{
		/* already in */
		return false;
	}

	AtomPtr c = seed->getAtom();
	if (std::find(atoms.begin(), atoms.end(), c) != atoms.end())
	{
		return true;
	}

	for (int i = 0; i < seed->connectedCount(); i++)
	{
		AtomPtr c = seed->getConnectedAtom(i);
		if (std::find(atoms.begin(), atoms.end(), c) != atoms.end())
		{
			/* connection somewhere */
			return true;
		}
	}
	
	/* no connection to list */
	return false;
}

bool addSpongeToList(SpongePtr seed, std::vector<AtomPtr> &atoms,
                     std::vector<SpongePtr> &sponges)
{
	bool added = false;

	if (std::find(sponges.begin(), sponges.end(), seed) == sponges.end())
	{
		sponges.push_back(seed);
	}

	AtomPtr c = seed->getAtom();
	if (std::find(atoms.begin(), atoms.end(), c) == atoms.end())
	{
		atoms.push_back(c);
	}

	for (int i = 0; i < seed->connectedCount(); i++)
	{
		AtomPtr c = seed->getConnectedAtom(i);
		if (std::find(atoms.begin(), atoms.end(), c) == atoms.end())
		{
			atoms.push_back(c);
		}
	}
	
	return added;
}

std::vector<SpongePtr> WaterNetwork::acquireSponges(SpongePtr seed)
{
	std::vector<AtomPtr> atoms;
	std::vector<SpongePtr> sponges;
	
	addSpongeToList(seed, atoms, sponges);
	
	bool changed = true;

	while (changed)
	{
		changed = false;
		for (int i = 0; i < atomCount(); i++)
		{
			if (atom(i)->getModel()->isSponge())
			{
				SpongePtr sponge = ToSpongePtr(atom(i)->getModel());
				
				if (spongeToBeAdded(sponge, atoms, sponges))
				{
					changed |= addSpongeToList(sponge, atoms, sponges);
					continue;
				}
			}
		}
	}

	return sponges;
}

void WaterNetwork::generateCalculations()
{
	if (_sponges.size() == 0)
	{
		return;
	}
	
	size_t num = _sponges[0]->sampleSize();

	_calcs = std::vector<WaterCalcs>(num, WaterCalcs());
	for (int i = 0; i < _sponges.size(); i++)
	{
		SpongePtr sponge = _sponges[i];
		
		for (int j = 0; j < num; j++)
		{
			WaterCalc calc;
			calc.middle = sponge->samplePointer(j);

			for (int k = 0; k < sponge->connectedCount(); k++)
			{
				AtomPtr con = sponge->getConnectedAtom(k);
				calc.left = con->getExplicitModel()->samplePointer(j);
				calc.ratio = sponge->distance(k);
				calc.weight = 0.1;
				
				if (!con->getModel()->isSponge())
				{
					calc.weight = 1;
				}

				_calcs[j].push_back(calc);
			}
		}
	}
}

void WaterNetwork::setActive(int a)
{
	for (int i = 0; i < _sponges.size(); i++)
	{
		_sponges[i]->setActive(a);
	}
}

double WaterNetwork::score()
{
	for (int i = 0; i < _sponges.size(); i++)
	{
		/* for display in GUI */
		_sponges[i]->copyActiveToFinalPos();
	}
	
	double sum = 0;

	for (int i = 0; i < _calcs[_n].size(); i++)
	{
		WaterCalc calc = _calcs[_n][i];
		vec3 diff = vec3_subtract_vec3(*(calc.left), *(calc.middle));
		double length = vec3_length(diff);
		double add = (length - calc.ratio);
		add *= add;
		add *= calc.weight;
		sum += add;
	}

	return sum;
}

void WaterNetwork::calculateSingle(SpongePtr sponge, bool others)
{
	_sponges.clear();
	_sponges.push_back(sponge);
	
	if (others)
	{
		_sponges = acquireSponges(sponge);
	}
	
	generateCalculations();
	refineSponges();
}

void WaterNetwork::refineSponges()
{
	double step = 0.01;
	double tol = 0.0005;
	bool changed = true;
	int num = 0;
	
	while (changed && num < 100)
	{
		num++;
		changed = false;

		for (int i = 0; i < _calcs.size(); i++)
		{
			_n = i;
			setActive(i);

			NelderMeadPtr neld = NelderMeadPtr(new RefinementNelderMead());
			neld->setJobName("sponge");
			neld->setSilent();
			neld->setCycles(200);

			for (int j = 0; j < _sponges.size(); j++)
			{
				neld->addParameter(&*_sponges[j], Sponge::getX, Sponge::setX,
				                   step, tol, "s" + i_to_str(i) + "_x");
				neld->addParameter(&*_sponges[j], Sponge::getY, Sponge::setY,
				                   step, tol, "s" + i_to_str(i) + "_y");
				neld->addParameter(&*_sponges[j], Sponge::getZ, Sponge::setZ,
				                   step, tol, "s" + i_to_str(i) + "_z");
			}

			neld->setEvaluationFunction(WaterNetwork::sScore, this);
			neld->refine();
			changed |= neld->changedSignificantly();
		}
	}

	for (int i = 0; i < _sponges.size(); i++)
	{
		_sponges[i]->recalculated();
	}

}

SpongePtr WaterNetwork::findFirstSponge()
{
	for (int i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getModel()->isSponge())
		{
			SpongePtr sponge = ToSpongePtr(atom(i)->getModel());
			if (sponge->needsRecalc())
			{
				return sponge;
			}
		}
	}
	
	return SpongePtr();
}

void WaterNetwork::refreshSponges()
{
	for (int i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getModel()->isSponge())
		{
			SpongePtr sponge = ToSpongePtr(atom(i)->getModel());
			sponge->propagateChange();
			sponge->Novalent::getManyPositions();
			sponge->propagateChange();
			calculateSingle(sponge, false);
			sponge->propagateChange();
		}
	}
}

void WaterNetwork::recalculate()
{
	refreshSponges();
	SpongePtr sponge = findFirstSponge();

	while (sponge)
	{
		_sponges = acquireSponges(sponge);
		std::cout << "Acquired " << _sponges.size() << " sponges" << std::endl;
		
		if (_sponges.size())
		{
			generateCalculations();
			refineSponges();
			sponge = findFirstSponge();
		}
		else
		{
			break;
		}
	}
}

void WaterNetwork::prune()
{
	int deleted = 0;
	int total = atomCount();
	std::cout << "Pruning waters" << std::endl;
	for (int i = 0; i < atomCount(); i++)
	{
		if (!atom(i)->getModel()->isAbsolute())
		{
			continue;
		}

		MapScoreWorkspace ws;
		setup_space(&ws);
		ws.crystal = Options::getActiveCrystal();

		AtomGroupPtr single = AtomGroupPtr(new AtomGroup());
		single->addAtom(atom(i));
		AtomList extra = ws.crystal->getCloseAtoms(atom(i), 3.0, false);
		single->addAtomsFrom(extra);
		
		ws.selectAtoms = single;
		
		double score = AtomGroup::scoreWithMapGeneral(&ws);
		double orig = atom(i)->getModel()->getEffectiveOccupancy();
		ToAbsolutePtr(atom(i)->getModel())->setOccupancy(0.1);
		double newscore = AtomGroup::scoreWithMapGeneral(&ws);
		
		std::cout << score << " to " << newscore << std::flush;
		if (newscore < score)
		{
			std::cout << "- deleted" << std::endl;
			ws.crystal->removeAtom(atom(i));
			deleted++;
			i--;
		}
		else
		{
			ToAbsolutePtr(atom(i)->getModel())->setOccupancy(orig);
			std::cout << std::endl;
		}
	}

	std::cout << "Deleted " << deleted << " out of " << total
	<< " waters (" << (double)deleted/(double)total *100 << "%)." 
	<< std::endl;
}


