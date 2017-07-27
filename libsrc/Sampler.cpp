//
//  Sampler.cpp
//  vagabond
//
//  Created by Helen Ginn on 25/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Sampler.h"
#include "RefinementGridSearch.h"
#include "Bond.h"
#include "Atom.h"
#include "Crystal.h"
#include "Element.h"

Sampler::Sampler()
{

}

void Sampler::addTorsion(BondPtr bond, double range, double interval)
{
	if (!_strategy)
	{
		_strategy = RefinementStrategyPtr(new RefinementGridSearch());
		_strategy->setEvaluationFunction(Sampler::score, this);
	}

	double number = fabs(range / interval);
	_strategy->addParameter(&*bond, Bond::getTorsion, Bond::setTorsion,
							interval, number, "torsion");

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
	_strategy->setJobName(_jobName);
	_strategy->refine();
	_strategy = RefinementStrategyPtr();

	double rad = Bond::getTorsion(&*_bonds[0]);
	std::cout << "Torsion: " << rad2deg(rad) << std::endl;
}


double Sampler::getScore()
{
	double score = 1;

	for (int i = 0; i < _sampled.size(); i++)
	{
		double next_score = _sampled[i]->scoreWithMap(_fft, _real2hkl);
		score *= next_score;
	}

	return -score;
}