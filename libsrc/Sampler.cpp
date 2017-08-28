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
#include "Polymer.h"
#include "Monomer.h"
#include "Backbone.h"

Sampler::Sampler()
{
	_mock = false;
	_joint = false;
	_scoreType = ScoreTypeCorrel;
}


void Sampler::setupTorsionSet(BondPtr bond, int k, int bondNum, int resNum,
								 double range, double interval, bool addDampen)
{
	bond->setActiveGroup(k);
	bond->setBlocked(true);

	reportInDegrees();
	setScoreType(ScoreTypeRFactor);

	setJobName("torsion_set_" + bond->getMajor()->getAtomName() + "_" +
			   bond->getMinor()->getAtomName() + "_g" +
			   i_to_str(k) + "_" + i_to_str(resNum));

	if (addDampen)
	{
		addTorsionBlur(bond, 0.01, 0.1);
		addDampening(bond, 0.2, 0.01);
	}
	else
	{
		addTorsion(bond, deg2rad(range), deg2rad(interval));
	}

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
				if (!addDampen)
				{
					addTorsion(nextBond, deg2rad(range), deg2rad(interval));
				}
				else
				{
					addDampening(nextBond, 0.2, 0.01);
				}

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

	if (addDampen)
	{
		addDampening(bond, 0.2, 0.01);
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
	_strategy->setCycles(20);

}

void Sampler::addOccupancy(BondPtr bond, double range, double interval)
{
	//	double number = fabs(range / interval);
	_strategy->addParameter(&*bond, Bond::getOccupancy, Bond::setOccupancy,
							range, interval,
							"occupancy");

	_bonds.push_back(bond);
}

void Sampler::addRamachandranAngles(PolymerPtr polymer, int from, int to)
{
	for (int i = from - 1; i < to - 1; i++)
	{
		if (!polymer->getMonomer(i))
		{
			continue;
		}

		BackbonePtr backbone = polymer->getMonomer(i)->getBackbone();
		AtomPtr ca = backbone->findAtom("CA");
		AtomPtr n = backbone->findAtom("C");

		if (ca)
		{
			if (ca->getModel()->getClassName() != "Bond")
			{
				continue;
			}

			BondPtr caBond = ToBondPtr(ca->getModel());
			addTorsion(caBond, deg2rad(0.2), deg2rad(0.05));

			if (caBond->getParentModel()->getClassName() == "Absolute")
			{
				AbsolutePtr abs = ToAbsolutePtr(caBond->getParentModel());
				addAbsolutePosition(abs, 0.02, 0.01);
			}


		}

		if (n)
		{
			if (n->getModel()->getClassName() != "Bond")
			{
				continue;
			}

			BondPtr nBond = ToBondPtr(n->getModel());
			addTorsion(nBond, deg2rad(0.2), deg2rad(0.05));
		}
	}
}

void Sampler::addTorsion(BondPtr bond, double range, double interval)
{
//	double number = fabs(range / interval);
	std::string num = i_to_str(_strategy->parameterCount() + 1);
	_strategy->addParameter(&*bond, Bond::getTorsion, Bond::setTorsion,
							range, interval,
							"t" + bond->shortDesc());

	_bonds.push_back(bond);
}

void Sampler::addMagicAxis(BondPtr bond, double range, double interval)
{
	//	double number = fabs(range / interval);
	std::string num = i_to_str(_strategy->parameterCount() + 1);

	_strategy->addParameter(&*bond, Bond::getHRot, Bond::setHRot,
							range, interval,
							"h" + bond->shortDesc());
	_strategy->addParameter(&*bond, Bond::getKRot, Bond::setKRot,
							range, interval,
							"k" + bond->shortDesc());

//	_strategy->setSilent(true);
	_bonds.push_back(bond);
}

void Sampler::addMagicAxisBroad(BondPtr bond)
{
	_strategy->addParameter(&*bond, Bond::getMagicAxisMat, Bond::setMagicAxisMat,
							6.0, 1.01,
							"ax" + bond->shortDesc());
	_bonds.push_back(bond);
}

void Sampler::addTorsionBlur(BondPtr bond, double range, double interval)
{
	_strategy->addParameter(&*bond, Bond::getTorsionBlur, Bond::setTorsionBlur,
							range, interval, "b" + bond->shortDesc());
	_strategy->addParameter(&*bond, Bond::getTorsionVertBlur, Bond::setTorsionVertBlur,
							range, interval, "vb" + bond->shortDesc());

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

void Sampler::addDampening(BondPtr bond, double range, double interval)
{
//	double number = fabs(range / interval);
	_strategy->addParameter(&*bond, Bond::getDampening,
							Bond::setDampening, range,
							interval, "d" + bond->shortDesc());

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

void Sampler::addSampledCAs(PolymerPtr polymer, int from, int to)
{
	for (int i = from; i < to; i++)
	{
		if (!polymer->getMonomer(i))
		{
			continue;
		}

		BackbonePtr backbone = polymer->getMonomer(i)->getBackbone();
		AtomPtr ca = backbone->findAtom("CA");
		addSampled(ca);
		AtomPtr c = backbone->findAtom("C");
		addSampled(c);
		AtomPtr n = backbone->findAtom("N");
		addSampled(n);
	}
}

void Sampler::addSampledAtoms(AtomGroupPtr group)
{
	if (!group)
	{
		return;
	}

	for (int i = 0; i < group->atomCount(); i++)
	{
		addSampled(group->atom(i));
	}
}

void Sampler::addSampled(AtomPtr atom)
{
	if (!atom)
	{
		return;
	}

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

	if (_scoreType == ScoreTypeModelPos)
	{
		_strategy->setCycles(50);
	}

	if (sampleSize())
	{
		_strategy->setJobName(_jobName);
		_strategy->refine();
	}

	_scoreType = ScoreTypeCorrel;
	_strategy = RefinementStrategyPtr();

	for (int i = 0; i < _bonds.size(); i++)
	{
		_bonds[i]->setBlocked(false);
	}

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

	if (_scoreType == ScoreTypeModelRMSD || _scoreType == ScoreTypeModelRMSDZero)
	{
		double score = 0;
		double count = 0;

		for (int i = 0; i < sampleSize(); i++)
		{
			ModelPtr model = _sampled[i]->getModel();

			if (model->getClassName() != "Bond")
			{
				continue;
			}

			BondPtr bond = std::static_pointer_cast<Bond>(model);
			double target = -1;

			if (_scoreType == ScoreTypeModelRMSD)
			{
				target = _sampled[i]->getInitialBFactor();
			}

			double rmsdScore = bond->getMeanSquareDeviation(target);

			count++;
			score += rmsdScore;
		}

		return score / count;
	}

	if (_scoreType == ScoreTypeModelPos)
	{
		double score = 0;
		
		for (int i = 0; i < sampleSize(); i++)
		{
			BondPtr bond = ToBondPtr(_sampled[i]->getModel());
			bond->getDistribution();
			vec3 bestPos = bond->getAbsolutePosition();
			vec3 initialPos = _sampled[i]->getInitialPosition();

			vec3 diff = vec3_subtract_vec3(bestPos, initialPos);
			score += vec3_sqlength(diff);
		}

		score /= (double)sampleSize();

		return score;
	}

	std::vector<double> xs, ys;

	double n = 60;
	double scales = 0.33;
	FFTPtr segment = FFTPtr(new FFT());
	segment->create(n);
	segment->setScales(scales);
	mat3x3 basis = make_mat3x3();
	double toReal = 1/(scales*n);
	mat3x3_scale(&basis, toReal, toReal, toReal);

	_sampled[0]->getModel()->getDistribution();
	vec3 offset = _sampled[0]->getModel()->getAbsolutePosition();

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
	else if (_scoreType == ScoreTypeRFactor)
	{
		double rFactor = scaled_r_factor(xs, ys, cutoff);
		return rFactor;
	}
	else
	{
		double mult = weightedMapScore(xs, ys);
		return -mult;
	}
}