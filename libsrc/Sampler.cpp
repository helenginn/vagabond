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
#include "Absolute.h"

Sampler::Sampler()
{
	_mock = false;
	_joint = false;
	_scoreType = ScoreTypeCorrel;
}


void Sampler::setupDoubleTorsion(BondPtr bond, int k, int bondNum, int resNum,
								 double range, double interval)
{
	bond->setActiveGroup(k);

	setupGrid();
	reportInDegrees();

	setJobName("torsion_double_" + bond->getMajor()->getAtomName() + "_" +
			   bond->getMinor()->getAtomName() + "_g" +
			   i_to_str(k) + "_" + i_to_str(resNum));

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

		int trial = 0;
		AtomPtr nextAtom = bond->downstreamAtom(k, trial);
		int totalAtoms = bond->downstreamAtomCount(k);

		BondPtr nextBond = std::static_pointer_cast<Bond>(nextAtom->getModel());

		while (nextAtom &&
			   (!nextBond->isNotJustForHydrogens() || nextBond->isFixed() ||
				!nextBond->isUsingTorsion()))
		{
			trial++;
			if (trial >= totalAtoms)
			{
				break;
			}
			nextAtom = bond->downstreamAtom(k, trial);
			nextBond = std::static_pointer_cast<Bond>(nextAtom->getModel());
		}

		if (nextAtom)
		{
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
							bond->shortDesc() + "-torsion");

	_bonds.push_back(bond);
}

void Sampler::addTorsionBlur(BondPtr bond, double range, double interval)
{
	_strategy->addParameter(&*bond, Bond::getTorsionBlur, Bond::setTorsionBlur,
							range, interval, bond->shortDesc() + "-tblur");

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
	//	double number = fabs(range / interval);
	_strategy->addParameter(&*bond, Bond::getBendAngle, Bond::setBendAngle,
							range, interval, "bend");

	_bonds.push_back(bond);
}


void Sampler::addAbsolutePosition(AbsolutePtr abs, double range, double interval)
{
	if (abs->getClassName() != "Absolute")
	{
		return;
	}

	//	double number = fabs(range / interval);
	_strategy->addParameter(&*abs, Absolute::getPosX, Absolute::setPosX,
							range, interval, "pos_x");
	_strategy->addParameter(&*abs, Absolute::getPosY, Absolute::setPosY,
							range, interval, "pos_y");
	_strategy->addParameter(&*abs, Absolute::getPosZ, Absolute::setPosZ,
							range, interval, "pos_z");
}


void Sampler::addAbsoluteBFactor(AbsolutePtr abs, double range, double interval)
{
	if (abs->getClassName() != "Absolute")
	{
		return;
	}

	//	double number = fabs(range / interval);
	_strategy->addParameter(&*abs, Absolute::getB, Absolute::setB,
							range, interval, "bfactor");
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
		for (int i = 0; i < sampleSize(); i++)
		{
			if (_sampled[i] == atom)
			{
				return;
			}
		}
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
/*
	std::cout << "Refining bonds ";

	for (int i = 0; i < _bonds.size(); i++)
	{
		std::cout << _bonds[i]->shortDesc() << " (" <<
		rad2deg(_bonds[i]->getTorsion(0)) << "º) ";
	}

	std::cout << "against atoms ";

	for (int i = 0; i < sampleSize(); i++)
	{
		std::cout << _sampled[i]->getAtomName() << " ";
	}

	std::cout << std::endl;*/

	if (sampleSize())
	{
		_strategy->setJobName(_jobName);
		_strategy->refine();
	}

	_scoreType = ScoreTypeCorrel;
	_strategy = RefinementStrategyPtr();
	_bonds.clear();
	_sampled.clear();
	_unsampled.clear();
	_joint = false;
}


double Sampler::getScore()
{
	if (!_sampled.size())
	{
		return 0;
	}

	std::vector<double> xs, ys;

	double n = 32;
	double scales = 0.33;
	FFTPtr segment = FFTPtr(new FFT());
	segment->create(n);
	segment->setScales(scales);
	mat3x3 basis = make_mat3x3();
	double toReal = 1/(scales*n);
	mat3x3_scale(&basis, toReal, toReal, toReal);

	vec3 offset = _sampled[0]->getPosition();

	for (int i = 0; i < _sampled.size(); i++)
	{
		_sampled[i]->addToMap(segment, basis, offset);
	}

	mat3x3_mult_vec(_real2hkl, &offset);
	double cutoff = FFT::score(_fft, segment, offset, &xs, &ys);

	if (_scoreType == ScoreTypeCorrel)
	{
		double correl = correlation(xs, ys, cutoff);
		return -correl;
	}
	else
	{
		double mult = weightedMapScore(xs, ys);
		return -mult;
	}
}