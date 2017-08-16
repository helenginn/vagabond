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
#include "FileReader.h"

Sampler::Sampler()
{
	_mock = false;
	_joint = false;
}


void Sampler::setupDoubleTorsion(BondPtr bond, int k, int bondNum, int resNum,
								 double range, double interval)
{
	bond->setActiveGroup(k);

	setupGrid();
	reportInDegrees();

	addTorsion(bond, deg2rad(range), deg2rad(interval));

	for (int j = 0; j < bond->downstreamAtomCount(k); j++)
	{
		addSampled(bond->downstreamAtom(k, j));
	}

	for (int j = 0; j < bond->extraTorsionSampleCount(k); j++)
	{
		addSampled(bond->extraTorsionSample(k, j));
	}

	for (int i = 0; i < bondNum; i++)
	{
		if (!bond->downstreamAtomGroupCount() || !bond->downstreamAtomCount(0) ||
			!bond->isUsingTorsion() || bond->isFixed() || !bond->isNotJustForHydrogens())
		{
			break;
		}

		AtomPtr nextAtom = bond->downstreamAtom(k, 0);
		BondPtr nextBond;

		if (nextAtom)
		{
			nextBond = std::static_pointer_cast<Bond>(nextAtom->getModel());

			if (nextBond->isUsingTorsion() && !nextBond->isFixed()
				&& nextBond->isNotJustForHydrogens())
			{
				addTorsion(nextBond, deg2rad(range), deg2rad(interval));

				for (int j = 0; j < nextBond->downstreamAtomCount(0); j++)
				{
					addSampled(nextBond->downstreamAtom(0, j));
				}

				for (int j = 0; j < nextBond->extraTorsionSampleCount(0); j++)
				{
					addSampled(nextBond->extraTorsionSample(0, j));
				}
			}
		}

		if (!nextBond)
		{
			break;
		}

		bond = nextBond;
	}

	setJobName("torsion_double_" + bond->getMajor()->getAtomName() + "_" +
			   bond->getMinor()->getAtomName() + "_g" +
			   i_to_str(k) + "_" + i_to_str(resNum));
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

void Sampler::addOccupancy(BondPtr bond, double range, double interval)
{
	//	double number = fabs(range / interval);
	_strategy->addParameter(&*bond, Bond::getOccupancy, Bond::setOccupancy,
							range, interval,
							"occupancy");

	_bonds.push_back(bond);
}

void Sampler::addTorsion(BondPtr bond, double range, double interval)
{
//	double number = fabs(range / interval);
	std::string num = i_to_str(_strategy->parameterCount() + 1);
	_strategy->addParameter(&*bond, Bond::getTorsion, Bond::setTorsion,
							range, interval,
							"torsion_" + num);

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
	_bonds.clear();
	_sampled.clear();
	_unsampled.clear();
	_joint = false;
}


double Sampler::getScore()
{
	double score = 1;

	std::vector<double> xTots, yTots;
	std::vector<double> xs, ys;

	std::vector<double> *xPtr = NULL;
	std::vector<double> *yPtr = NULL;

	if (_joint)
	{
		xPtr = &xs;
		yPtr = &ys;
	}

	for (int i = 0; i < _sampled.size(); i++)
	{
		double next_score = _sampled[i]->scoreWithMap(_fft, _real2hkl, xPtr, yPtr);

		next_score += 1;

		score *= next_score;

		if (_joint)
		{
			xTots.reserve(xTots.size() + xs.size());
			yTots.reserve(yTots.size() + ys.size());
			xTots.insert(xTots.end(), xs.begin(), xs.end());
			yTots.insert(yTots.end(), ys.begin(), ys.end());
		}
	}

	if (_joint)
	{
		double val = correlation(xTots, yTots);
		return -val;
	}

	return -fabs(score);
}