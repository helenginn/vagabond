//
//  Sampler.cpp
//  vagabond
//
//  Created by Helen Ginn on 25/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Sampler.h"
#include "RefinementGridSearch.h"
#include "RefinementNelderMead.h"
#include "Bond.h"
#include "Atom.h"
#include "Crystal.h"
#include "Element.h"

Sampler::Sampler()
{
	_mock = false;
}

void Sampler::setupGrid()
{
	_strategy = RefinementStrategyPtr(new RefinementGridSearch());
	_strategy->setEvaluationFunction(Sampler::score, this);
}

void Sampler::setupNelderMead()
{
	_strategy = RefinementStrategyPtr(new NelderMead());
	_strategy->setEvaluationFunction(Sampler::score, this);
	_strategy->setCycles(10);

}

void Sampler::addTorsion(BondPtr bond, double range, double interval)
{
//	double number = fabs(range / interval);
	_strategy->addParameter(&*bond, Bond::getTorsion, Bond::setTorsion,
							range, interval,
							"torsion");

	_bonds.push_back(bond);
}

void Sampler::addTorsionBlur(BondPtr bond, double range, double interval)
{
//	double number = fabs(range / interval);
	_strategy->addParameter(&*bond, Bond::getTorsionBlur, Bond::setTorsionBlur,
							range, interval, "torsion_blur");

	double blurNow = Bond::getTorsionBlur(&*bond);

	if (blurNow < range / 2)
	{
		Bond::setTorsionBlur(&*bond, range / 2);
	}

	_bonds.push_back(bond);
}

void Sampler::addBondLength(BondPtr bond, double range, double interval)
{
	//	double number = fabs(range / interval);
	_strategy->addParameter(&*bond, Bond::getBondLength,
							Bond::setBondLength, range,
							interval, "bond_length");

	_bonds.push_back(bond);
}

void Sampler::addTorsionNextBlur(BondPtr bond, double range, double interval)
{
//	double number = fabs(range / interval);
	_strategy->addParameter(&*bond, Bond::getTorsionNextBlur,
							Bond::setTorsionNextBlur, range,
							interval, "torsion_next_blur");

	_bonds.push_back(bond);
}

void Sampler::addBendBlur(BondPtr bond, double range, double interval)
{
//	double number = fabs(range / interval);
	_strategy->addParameter(&*bond, Bond::getBendBlur, Bond::setBendBlur,
							range, interval, "bend_blur");

	double blurNow = Bond::getBendBlur(&*bond);

	if (blurNow < (range / 2))
	{
		Bond::setBendBlur(&*bond, range / 2);
	}

	_bonds.push_back(bond);
}

void Sampler::addBendAngle(BondPtr bond, double range, double interval)
{
	if (bond->getParentModel()->getClassName() != "Bond")
	{
		return;
	}

	//	double number = fabs(range / interval);
	_strategy->addParameter(&*bond, Bond::getBendAngle, Bond::setBendAngle,
							range, interval, "bend");

	_bonds.push_back(bond);
}

void Sampler::addSampled(AtomPtr atom)
{
	double electrons = atom->getElement()->electronCount();

	if (electrons <= 1)
	{
		// hydrogen
		_unsampled.push_back(atom);
	}
	else
	{
		_sampled.push_back(atom);
	}
}

void Sampler::setCrystal(CrystalPtr crystal)
{
	_real2hkl = crystal->getReal2Frac();
	_fft = crystal->getFFT();
}

void Sampler::sample()
{
	if (_mock)
	{
		_strategy->isMock();
		_mock = false;
	}
	_strategy->setJobName(_jobName);
	_strategy->refine();
	_strategy = RefinementStrategyPtr();
	_sampled.clear();
	_unsampled.clear();
}


double Sampler::getScore()
{
	double score = 1;

	for (int i = 0; i < _sampled.size(); i++)
	{
		double next_score = _sampled[i]->scoreWithMap(_fft, _real2hkl);
		score *= next_score;
	}

//	score = pow(score, 1 / (double)_sampled.size());

	return -score;
}