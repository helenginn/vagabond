//  WaterCluster.cpp
//  vagabond
//
//  Created by Helen Ginn on 01/07/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "WaterCluster.h"
#include "RefinementNelderMead.h"
#include "RefinementStepSearch.h"
#include "Crystal.h"
#include "Atom.h"
#include "Element.h"
#include "Absolute.h"
#include "Options.h"

WaterCluster::WaterCluster()
{
	_modifySample = -1;
}

WaterCluster::WaterCluster(WaterCluster &other)
{
	_modifySample = -1;
	_waters = other._waters;
	_pairs = other._pairs;
	_atoms = other._atoms;
}

ChromosomalPtr WaterCluster::makeCopy()
{
	ChromosomalPtr copy = ChromosomalPtr(new WaterCluster(*this));
	
	return copy;
}

void WaterCluster::addAtom(AtomPtr atom)
{
	AtomGroup::addAtom(atom);

	if (atom->isHeteroAtom() && atom->getAtomName() == "O")
	{
		_waters.push_back(atom);
	}
}


void WaterCluster::findNeighbours()
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	int total = 0;
	int hBondy = 0;
	int allOpt = 0;

	for (int j = 0; j < _waters.size(); j++)
	{
		std::vector<AtomPtr> neighbours;
		neighbours = crystal->getCloseAtoms(_waters[j], 3.6);
		WaterNeighboursPair pair;
		pair.water = _waters[j];

		for (int i = 0; i < neighbours.size(); i++)
		{
			if (neighbours[i] == _waters[j])
			{
				continue;
			}
			
			Restraints options;
			options.push_back(RestraintNone);
			options.push_back(RestraintRepel);
			
			if (neighbours[i]->getElement()->getSymbol() == "H")
			{
				continue;
			}

			if (neighbours[i]->canBeHydrogenBonder())
			{
				hBondy++;
				options.push_back(RestraintHBond);
			}

			total++;
			RestraintChoice choice;
			choice.options = options;
			choice.picked = 1;
			RestraintMap map;
			map.neighbour = neighbours[i];
			map.choice = choice;
			double dist = _waters[j]->getDistanceFrom(&*neighbours[i], -1,
			                                          true);
			if (dist > 4.0) continue;

			map.close = true;
			
			pair.neighbours.push_back(map);

			allOpt += options.size();
		}

		_pairs.push_back(pair);
	}
	
	std::cout << "Cluster: " << _waters.size() << " waters have " 
	<< total << " neighbours. " << hBondy << " have "\
	"hydrogen bonding potential, and therefore " << allOpt << " binary"\
	" restraint options." << std::endl;
}

void WaterCluster::wipeBonding()
{
	for (int i = 0; i < _pairs.size(); i++)
	{
		for (int j = 0; j < _pairs[i].neighbours.size(); j++)
		{
			RestraintChoice *choice = &_pairs[i].neighbours[j].choice;
			choice->picked = 0;
		}
	}
}

double WaterCluster::scoreAgainstDensity()
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	double score;
	score = scoreWithMap(ScoreTypeCorrel, crystal);
	return score;
}

void WaterCluster::resetModels()
{
	for (int i = 0; i < _waters.size(); i++)
	{
		ModelPtr model = _waters[i]->getModel();
		if (model->isAbsolute())
		{
//			ToAbsolutePtr(model)->setImplicitPositions();
		}
	}
}

void WaterCluster::refine()
{
	findNeighbours();
	wipeBonding();
	double first = scoreAgainstDensity();
	std::cout << "B factor CC: " << first << std::endl;
	double second = evaluate();
	std::cout << "No bonding CC: " << second << std::endl;
	evolve();
	recalculateWaters();
	double third = scoreAgainstDensity();

	if (third < second)
	{
		std::cout << "Keeping improvement." << std::endl;
	}
	else if (second > first)
	{
		std::cout << "Resetting entirely." << std::endl;
		resetModels();
	}
	else
	{
		std::cout << "Keeping explicit positions only." << std::endl;
	}

	evaluate();
}

double WaterCluster::recalculateWaters()
{
	/*
	// Make sure all waters are an ensemble //
	for (int i = 0; i < _waters.size(); i++)
	{
		ModelPtr model = _waters[i]->getModel();
		if (!model->isAbsolute()) continue;

		AbsolutePtr abs = ToAbsolutePtr(model);
		
		if (abs->hasExplicitPositions())
		{
			// recalculate a fresh dist //
			abs->resetSamples();
			continue;
		}
		
		FFTPtr fft = abs->getDistribution();
		abs->setAnchorPoint();
		abs->refreshPositions();
	}
	
	int nSamples = Options::getRuntimeOptions()->getNSamples();
	
	for (int i = 0; i < nSamples; i++)
	{
		_modifySample = i;
		RefinementStepSearchPtr nelderMead = RefinementStepSearchPtr(new RefinementStepSearch());
		nelderMead->setCycles(100);
		nelderMead->setSilent(true);
		nelderMead->setEvaluationFunction(WaterCluster::score, this);

		// Make sure we only affect one of our samples //
		for (int j = 0; j < _waters.size(); j++)
		{
			ModelPtr model = _waters[j]->getModel();
			if (!model->isAbsolute()) continue;

			AbsolutePtr abs = ToAbsolutePtr(model);
			abs->setModifiedSample(i);	

			nelderMead->addParameter(&*abs, Absolute::getPosX,
			                         Absolute::setPosX,
			                         0.02, 0.002, "w" + i_to_str(j) + "x");
			nelderMead->addParameter(&*abs, Absolute::getPosY,
			                         Absolute::setPosY,
			                         0.02, 0.002, "w" + i_to_str(j) + "y");
			nelderMead->addParameter(&*abs, Absolute::getPosZ,
			                         Absolute::setPosZ,
			                         0.02, 0.002, "w" + i_to_str(j) + "z");
		}

		nelderMead->refine();
	}

	for (int j = 0; j < _waters.size(); j++)
	{
		ModelPtr model = _waters[j]->getModel();
		if (!model->isAbsolute()) continue;

		AbsolutePtr abs = ToAbsolutePtr(model);
		abs->clearModifiedSample();	
	}
	*/
}

double WaterCluster::evaluateRestraint(int sample, int i)
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	double contrib = 0;

	AtomPtr water = _pairs[i].water;

	std::vector<AtomPtr> allNeighs;

	for (int j = 0; j < _pairs[i].neighbours.size(); j++)
	{
		RestraintMap map = _pairs[i].neighbours[j];
		RestraintChoice choice = map.choice;

		RestraintType type = choice.options[choice.picked];
		if (type == RestraintNone)
		{
			continue;
		}
		
		AtomPtr neighbour = map.neighbour;
		bool quick = map.close;
		
		if (type == RestraintRepel)
		{
			double dist = water->getDistanceFrom(&*neighbour, sample, quick);
			double addition = 1 / dist;
			contrib += addition;
		}
		else if (type == RestraintHBond)
		{
			allNeighs.push_back(neighbour);
			double dist = water->getDistanceFrom(&*neighbour, sample, quick);
			double diff = 2.8 - dist;
			diff *= diff;
			contrib += diff;
		}
	}
	
	if (allNeighs.size() <= 1)
	{
		return contrib;		
	}
	
	// Hydrogen bond angles */
	/*
	vec3 waterPos = water->getModel()->getSpecificPosition(sample);
	double expected = deg2rad(104.5);

	for (int i = 1; i < allNeighs.size(); i++)
	{
		AtomPtr iN = allNeighs[i];
		vec3 iPos = iN->getModel()->getSpecificPosition(sample);
		iPos = vec3_subtract_vec3(iPos, waterPos);

		for (int j = 0; j < i; j++)
		{
			AtomPtr jN = allNeighs[j];
			vec3 jPos = jN->getModel()->getSpecificPosition(sample);
			jPos = vec3_subtract_vec3(jPos, waterPos);
			
			double angle = vec3_angle_with_vec3(iPos, jPos);
			double diff = angle - expected;
			diff *= 10;
			diff *= diff;
			contrib += diff;
		}
	}
	*/
	
	return contrib;
}

double WaterCluster::evaluateRestraints()
{
	double total = 0;
	
	for (int i = 0; i < _pairs.size(); i++)
	{
		double contrib = evaluateRestraint(_modifySample, i);
		total += contrib;
	}
	
	return total;
}

void WaterCluster::mutate()
{
	/* mutate 20% of choices */
	randomise(0.2);
}

void WaterCluster::randomise(double frac)
{
	for (int i = 0; i < _pairs.size(); i++)
	{
		for (int j = 0; j < _pairs[i].neighbours.size(); j++)
		{
			if (frac < 1 && rand() / (double)RAND_MAX > frac)
			{
				continue;
			}
			
			RestraintChoice *choice = &_pairs[i].neighbours[j].choice;
			size_t num = choice->options.size();
			int new_value = rand() % num;
			new_value *= (rand() % 2);

			choice->picked = new_value;
		}
	}
}

double WaterCluster::evaluate()
{
	recalculateWaters();
	double evaluation = scoreAgainstDensity();
	return evaluation;
}

void WaterCluster::haveSexWith(Chromosomal *_other)
{
	/* Swap over at some point in the waters */
	int pairCount = _pairs.size();
	int swapPoint = rand() % pairCount;
	WaterCluster *other = static_cast<WaterCluster *>(_other);

	for (int i = 0; i < swapPoint; i++)
	{
		for (int j = 0; j < _pairs[i].neighbours.size(); j++)
		{
			RestraintChoice *choice = &_pairs[i].neighbours[j].choice;
			int new_value = other->_pairs[i].neighbours[j].choice.picked;
			choice->picked = new_value;
		}
	}

	for (int i = swapPoint; i < _pairs.size(); i++)
	{
		for (int j = 0; j < _pairs[i].neighbours.size(); j++)
		{
			RestraintChoice *choice = &other->_pairs[i].neighbours[j].choice;
			int new_value = _pairs[i].neighbours[j].choice.picked;
			choice->picked = new_value;
		}
	}
}

void WaterCluster::geneticCode()
{
	for (int i = 0; i < _pairs.size(); i++)
	{
		for (int j = 0; j < _pairs[i].neighbours.size(); j++)
		{
			RestraintChoice *choice = &_pairs[i].neighbours[j].choice;
			std::cout << choice->picked << std::flush;
		}

		std::cout << "-" << std::flush;
	}

	std::cout << std::endl;
}

void WaterCluster::copyOver(Chromosomal *_other)
{
	WaterCluster *other = static_cast<WaterCluster *>(_other);

	for (int i = 0; i < _pairs.size(); i++)
	{
		for (int j = 0; j < _pairs[i].neighbours.size(); j++)
		{
			RestraintChoice *choice = &_pairs[i].neighbours[j].choice;
			int new_value = other->_pairs[i].neighbours[j].choice.picked;
			choice->picked = new_value;
		}
	}
}
