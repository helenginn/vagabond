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
#include "RefinementSnake.h"
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
	_refinedMagicAxisCount = 0;
	_overallFlex = 0.025;
}


BondPtr Sampler::setupTorsionSet(BondPtr bond, int k, int bondNum,
								 double range, double interval, bool addDampen)
{
	bond->setActiveGroup(k);

	reportInDegrees();
	setScoreType(ScoreTypeCorrel);

	setJobName("torsion_set_" + bond->shortDesc() + "_g" +
			   i_to_str(k));

	addTorsion(bond, deg2rad(range), deg2rad(interval));

	for (int j = 0; j < bond->downstreamAtomCount(k); j++)
	{
		addSampled(bond->downstreamAtom(k, j));
	}

	for (int j = 0; j < bond->extraTorsionSampleCount(k); j++)
	{
		addSampled(bond->extraTorsionSample(k, j));
	}

	BondPtr returnBond = BondPtr();

	for (int i = 0; i < bondNum; i++)
	{
		if (!bond->downstreamAtomGroupCount() || !bond->downstreamAtomCount(k) ||
			!bond->isRefinable())
		{
			/* No hope! Give up! */
			break;
		}

		int trial = 0;

		AtomPtr nextAtom = bond->downstreamAtom(k, trial);
		int totalAtoms = bond->downstreamAtomCount(k);

		BondPtr nextBond = std::static_pointer_cast<Bond>(nextAtom->getModel());

		/* In case the downstream atom has no future, get whatever
		 downstream atom has a refinable future.
		 Prioritise 0th bond. */
		while (nextAtom && !nextBond->isRefinable())
		{
			trial++;

			/* No more, give up */
			if (trial >= totalAtoms)
			{
				break;
			}

			nextAtom = bond->downstreamAtom(k, trial);
			nextBond = std::static_pointer_cast<Bond>(nextAtom->getModel());
		}

		/* Now if we have a better bond, go forth */
		if (nextAtom && nextBond->isRefinable())
		{
			if (!returnBond)
			{
				returnBond = nextBond;
			}

			addTorsion(nextBond, deg2rad(range), deg2rad(interval));

			for (int j = 0; j < nextBond->downstreamAtomCount(k); j++)
			{
				addSampled(nextBond->downstreamAtom(k, j));
			}

			for (int j = 0; j < nextBond->extraTorsionSampleCount(k); j++)
			{
				addSampled(nextBond->extraTorsionSample(k, j));
			}
		}

		if (!nextBond->isRefinable())
		{
			break;
		}

		bond = nextBond;
	}

	return returnBond;
}

void Sampler::setupSnake()
{
	_strategy = RefinementStrategyPtr(new RefinementSnake());
	RefinementSnakePtr snake = std::static_pointer_cast<RefinementSnake>(_strategy);
	snake->setParentSampler(this);
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

void Sampler::addOverallKickAndDampen(PolymerPtr polymer)
{
	_strategy->addParameter(&*polymer, Polymer::getBackboneDampening, Polymer::setBackboneDampening, 0.05, 0.02, "dampen");
}

void Sampler::addSidechainDampen(PolymerPtr polymer)
{
//	_strategy->addParameter(&*polymer, Polymer::getSidechainDampening, Polymer::setSidechainDampening, 0.02, 0.01, "side_dampen");
//	_strategy->addParameter(&*polymer, Polymer::getSideKick, Polymer::setSideKick, 0.02, 0.01, "side_kick");
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
	int step = (from < to) ? 1 : -1;

	for (int i = from - 1; i != to - 1; i += step)
	{
		if (!polymer->getMonomer(i))
		{
			continue;
		}

		std::string ramaAtom = "N";
		std::string peptideAtom = "C";

		if (polymer->getAnchor() >= polymer->getMonomer(i)->getResidueNum())
		{
			ramaAtom = "C";
			peptideAtom = "N";
		}

		reportInDegrees();
		BackbonePtr backbone = polymer->getMonomer(i)->getBackbone();
		AtomPtr ca = backbone->findAtom("CA");
		AtomPtr rama = backbone->findAtom(ramaAtom);
		AtomPtr peptide = backbone->findAtom(peptideAtom);
		std::vector<AtomPtr> atoms;

		BondPtr caBond = ToBondPtr(ca->getModel());

		BondPtr ramaBond = ToBondPtr(rama->getModel());
		BondPtr peptideBond = ToBondPtr(peptide->getModel());


		bool last = (i == to - 1 - step);

		if (!last)
		{
			addTorsion(peptideBond, ANGLE_SAMPLING, deg2rad(0.05));
			addTorsion(ramaBond, ANGLE_SAMPLING, deg2rad(0.05));
			addTorsion(caBond, ANGLE_SAMPLING, deg2rad(0.05));
		}
		else
		{
			addTorsion(ramaBond, ANGLE_SAMPLING, deg2rad(0.05));
			addTorsion(caBond, ANGLE_SAMPLING, deg2rad(0.05));
		}
	}
}

void Sampler::addTorsion(BondPtr bond, double range, double interval)
{
	if (!bond || bond->getClassName() != "Bond")
	{
		return;
	}

	int resNum = bond->getMinor()->getMonomer()->getResidueNum();

//	double number = fabs(range / interval);
	std::string num = i_to_str(_strategy->parameterCount() + 1);
	_strategy->addParameter(&*bond, Bond::getTorsion, Bond::setTorsion,
							range, interval,
							"t" + bond->shortDesc());

	_bonds.push_back(bond);
}

void Sampler::addMagicAngle(BondPtr bond, double range, double interval)
{
	if (!bond) return;

	//	double number = fabs(range / interval);
	std::string num = i_to_str(_strategy->parameterCount() + 1);

	_strategy->addParameter(&*bond, Bond::getMagicAngle, Bond::setMagicAngle,
							range, interval,
							"h" + bond->shortDesc());

//	_strategy->setSilent(true);
	_bonds.push_back(bond);
}


void Sampler::addTorsionBlur(BondPtr bond, double range, double interval)
{
	if (!bond) return;

	_strategy->addParameter(&*bond, Bond::getTorsionBlur, Bond::setTorsionBlur,
							range, interval, "b" + bond->shortDesc());
	
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
	if (!bond) return;

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

void Sampler::addSampledSidechains(PolymerPtr polymer)
{
	for (int i = 0; i < polymer->monomerCount(); i++)
	{
		if (!polymer->getMonomer(i))
		{
			continue;
		}

		SidechainPtr sidechain = polymer->getMonomer(i)->getSidechain();
		addSampledAtoms(sidechain);
	}
}

void Sampler::addSampledBackbone(PolymerPtr polymer, int from, int to)
{
	if (from == 0 && to == 0)
	{
		from = 0;
		to = polymer->monomerCount();
	}

	int step = (from < to) ? 1 : -1;

	for (int i = from - 1; i != to - 1; i += step)
	{
		if (!polymer->getMonomer(i))
		{
			continue;
		}

		BackbonePtr backbone = polymer->getMonomer(i)->getBackbone();
		SidechainPtr sidechain = polymer->getMonomer(i)->getSidechain();

		AtomPtr ca = backbone->findAtom("CA");
		addSampled(ca);
		AtomPtr cb = sidechain->findAtom("CB");
		if (polymer->getMonomer(i)->getIdentifier() == "pro")
		{
			addSampledAtoms(sidechain);
		}
		else
		{
			addSampled(cb);
		}


		AtomPtr o = backbone->findAtom("O");
		addSampled(o);
		AtomPtr c = backbone->findAtom("C");
		addSampled(c);
		AtomPtr n = backbone->findAtom("N");
		addSampled(n);
	}
}

void Sampler::addSampled(std::vector<AtomPtr> atoms)
{
	for (int i = 0; i < atoms.size(); i++)
	{
		addSampled(atoms[i]);
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

bool Sampler::sample(bool clear)
{
	if (_mock)
	{
		_strategy->isMock();
		_mock = false;
	}

	if (_scoreType == ScoreTypeModelPos)
	{
		_strategy->setCycles(12);
	}
	else if (_scoreType == ScoreTypeModelRMSDZero)
	{
		_strategy->setCycles(10);
	}

	if (sampleSize())
	{
		_strategy->setJobName(_jobName);
		_strategy->refine();
	}

	_scoreType = ScoreTypeCorrel;

	bool changed = _strategy->didChange();

	if (clear)
	{
		_strategy = RefinementStrategyPtr();
	}

	for (int i = 0; i < _bonds.size(); i++)
	{
		_bonds[i]->setBlocked(false);
	}

	_bonds.clear();
	_sampled.clear();
	_unsampled.clear();
	_joint = false;

	return changed;
}

double Sampler::getScore()
{
	if (!_sampled.size())
	{
		return 0;
	}

	if (_scoreType == ScoreTypeModelFlexiness
		|| _scoreType == ScoreTypeModelRMSDZero)
	{
		double score = 0;
		double count = 0;

		for (int i = 0; i < sampleSize(); i++)
		{
			_sampled[i]->getModel()->propagateChange();
		}

		for (int i = 0; i < sampleSize(); i++)
		{
			ModelPtr model = _sampled[i]->getModel();

			if (model->getClassName() != "Bond")
			{
				continue;
			}

			BondPtr bond = std::static_pointer_cast<Bond>(model);
			double target = -1;

			if (_scoreType == ScoreTypeModelFlexiness)
			{
				target = _overallFlex;
			}

			double rmsdScore = bond->getFlexibilityPotential();

			count++;
			score += fabs(rmsdScore - target);
		}

		return score / count;
	}

	if (_scoreType == ScoreTypeModelPos)
	{
		double score = 0;
		double count = 0;

		for (int i = 0; i < sampleSize(); i++)
		{
			double oneScore = _sampled[i]->posDisplacement();
			score += oneScore;
			count++;
		}

		return score / count;

		for (int i = 1; i < sampleSize(); i++)
		{
			BondPtr bond = ToBondPtr(_sampled[i]->getModel());
			bond->getDistribution();
			vec3 bestPos = bond->getAbsolutePosition();
			vec3 initialPos = _sampled[i]->getPDBPosition();

			for (int j = 0; j < i - 1; j++)
			{
				BondPtr jBond = ToBondPtr(_sampled[j]->getModel());
				vec3 jBestPos = jBond->getAbsolutePosition();
				vec3 jInitialPos = _sampled[j]->getPDBPosition();

				// take (opposite) differences
				vec3 bestPosDiff = vec3_subtract_vec3(jBestPos, bestPos);
				vec3 initialDiff = vec3_subtract_vec3(initialPos, jInitialPos);

				// add them up and take the length of the remainder.
				vec3 addDiff = vec3_add_vec3(initialDiff, bestPosDiff);
				double length = vec3_length(addDiff);

				score += length;
				count++;
			}
		}

		score /= count;

		return score;
	}

	std::vector<double> xs, ys;

	_sampled[0]->getModel()->getDistribution();
	vec3 zero = _sampled[0]->getModel()->getAbsolutePosition();
	double maxDistance = 0;

	for (int i = 1; i < _sampled.size(); i++)
	{
		/* Refresh absolute position */
		_sampled[i]->getModel()->getDistribution();
		vec3 offset = _sampled[i]->getModel()->getAbsolutePosition();

		vec3 diff = vec3_subtract_vec3(offset, zero);
		double distance = vec3_length(diff);

		if (distance > maxDistance)
		{
			maxDistance = distance;
		}
	}

	double scales = 0.33;
	double n = (maxDistance + 1.0) / scales;
	n = 60;

	FFTPtr segment = FFTPtr(new FFT());
	segment->create(n + 0.5);
	segment->setScales(scales);
	mat3x3 basis = make_mat3x3();
	double toReal = 1/(scales*n);
	mat3x3_scale(&basis, toReal, toReal, toReal);

	for (int i = 0; i < _sampled.size(); i++)
	{
		_sampled[i]->addToMap(segment, basis, zero);
	}

	mat3x3_mult_vec(_real2hkl, &zero);
	double cutoff = FFT::score(_fft, segment, zero, &xs, &ys);

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

std::vector<double> Sampler::getNextResult(int num)
{
	RefinementGridSearchPtr grid;
	grid = std::static_pointer_cast<RefinementGridSearch>(_strategy);

	return grid->getNextResult(num);
}
