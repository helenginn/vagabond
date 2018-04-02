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

typedef struct
{
	BondPtr bond;
	int num;
} BondInt;

Sampler::Sampler()
{
	_mock = false;
	_scoreType = ScoreTypeCorrel;
	_refinedMagicAxisCount = 0;
	_silent = false;
}

void Sampler::addAtomsForBond(BondPtr firstBond, int k)
{
	std::vector<BondPtr> extraTorsionBonds;
	extraTorsionBonds.push_back(firstBond);

	for (int i = 0; i < extraTorsionBonds.size(); i++)
	{
		BondPtr bond = extraTorsionBonds[i];
		addSampled(bond->getMinor());

		for (int j = 0; j < bond->downstreamAtomCount(k); j++)
		{
			AtomPtr downAtom = bond->downstreamAtom(k, j);
			addSampled(downAtom);

			if (downAtom->getModel()->isBond())
			{
				BondPtr downBond = ToBondPtr(downAtom->getModel());
			}
		}

		for (int j = 0; j < bond->extraTorsionSampleCount(k); j++)
		{
			addSampled(bond->extraTorsionSample(k, j));
		}

		k = 0;
	}
}

void Sampler::addParamsForBond(BondPtr bond)
{
	for (ParamMap::iterator it = _params.begin(); it != _params.end(); it++)
	{
		ParamOptionType option = it->first;
		double range = it->second;

		switch (option)
		{
			case ParamOptionTorsion:
			addTorsion(bond, deg2rad(range), deg2rad(0.010));
			break;

			case ParamOptionKick:
			addTorsionBlur(bond, range, 0.01);
			break;

			case ParamOptionDampen:
			addDampening(bond, range, 0.0001);
			break;

			case ParamOptionMagicAngles:
			addMagicAngle(bond, deg2rad(range), deg2rad(0.10));
			break;

			default:
			break;
		}
	}
}

BondPtr Sampler::setupThoroughSet(BondPtr bond, int bondNum,
                                 double range, double interval, bool addAngle,
bool addFlex)
{
	reportInDegrees();
	setScoreType(ScoreTypeCorrel);

	setJobName("torsion_set_" + bond->shortDesc() + "_g");

	int bondCount = 0;

	if (_params.count(ParamOptionNumBonds))
	{
		bondNum = _params[ParamOptionNumBonds] + 0.2;
	}
	
	/* Logic: bond add map keeps remaining bonds (and remaining
 * 	additions) in memory for adding (allowing branches). When a
 * 	branched bond needs adding it is appended to this list with
 * 	one fewer remaining additions and the old entry deleted. */

	std::vector<BondInt> remaining;
	BondInt entry;
	entry.bond = bond;
	entry.num = bondNum;
	remaining.push_back(entry);
	addSampled(bond->getMajor());

	while (remaining.size())
	{
		BondInt first = remaining[0];
		BondPtr bond = first.bond;
		int num = first.num;
		
		remaining.erase(remaining.begin());

		if (num <= 0)
		{
			continue;	
		}
		
		bondCount++;
		addParamsForBond(bond);

		if (addAngle) 
		{
			addBendAngle(bond, deg2rad(0.01), deg2rad(0.001));
		}

		if (!bond->isRefinable())
		{
			/* No hope! Give up! */
			continue;
		}

		for (int j = 0; j < bond->downstreamAtomGroupCount(); j++)
		{
			addAtomsForBond(bond, j);

			/* Take the chosen group and check for futures */
			for (int i = 0; i < bond->downstreamAtomCount(j); i++)
			{
				AtomPtr downstreamAtom = bond->downstreamAtom(j, i);
				BondPtr nextBond = ToBondPtr(downstreamAtom->getModel());

				BondInt entry;
				entry.bond = nextBond;
				entry.num = num - 1;
				remaining.push_back(entry);
			}
		}
	}
	
	return BondPtr();
}
	
BondPtr Sampler::setupTorsionSet(BondPtr bond, int k, int bondNum,
                                 double range, double interval, bool
								 addAngle, bool addFlex)
{
	bond->setActiveGroup(k);

	reportInDegrees();
	setScoreType(ScoreTypeCorrel);

	setJobName("torsion_set_" + bond->shortDesc() + "_g" +
	           i_to_str(k));

	int bondCount = 0;

	if (_params.count(ParamOptionNumBonds))
	{
		bondNum = _params[ParamOptionNumBonds] + 0.2;
	}
	
	addSampled(bond->getMajor());

	if (addAngle) 
	{
		addBendAngle(bond, deg2rad(0.01), deg2rad(0.001));
	}

	addParamsForBond(bond);
	addAtomsForBond(bond, k);

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

		BondPtr nextBond = boost::static_pointer_cast<Bond>(nextAtom->getModel());

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
			nextBond = boost::static_pointer_cast<Bond>(nextAtom->getModel());
		}

		/* Now if we have a better bond, go forth */
		if (nextAtom && nextBond->isRefinable())
		{
			if (!returnBond)
			{
				returnBond = nextBond;
			}

			if (addAngle)
			{
				addBendAngle(nextBond, deg2rad(0.01), deg2rad(0.001));
			}

			addAtomsForBond(nextBond, 0);
			addParamsForBond(nextBond);
		}

		if (!nextBond->isRefinable())
		{
			break;
		}

		k = 0;
		bond = nextBond;
		bondCount++;
	}

	if (bondCount <= 2)
	{
		//    return BondPtr();
	}

	return returnBond;
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
	//    double number = fabs(range / interval);
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
	if (!bond || !bond->isBond())
	{
		return;
	}

	if (!bond->isRefinable())
	{
		return;
	}

	double mult = bond->getTorsionStepMult();
	range *= mult;

	_strategy->addParameter(&*bond, Bond::getTorsion, Bond::setTorsion,
	                        range, interval,
	"t" + bond->shortDesc());


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
	//    double number = fabs(range / interval);
	_strategy->addParameter(&*bond, Bond::getBondLength,
	                        Bond::setBondLength, range,
	interval, "bond_length");

	_bonds.push_back(bond);
}

void Sampler::addDampening(BondPtr bond, double range, double interval)
{
	if (!bond) return;

	//    double number = fabs(range / interval);
	_strategy->addParameter(&*bond, Bond::getDampening,
	                        Bond::setDampening, range,
	interval, "d" + bond->shortDesc());

	_bonds.push_back(bond);
}

void Sampler::addMagicAngle(BondPtr bond, double range, double interval)
{
	if (!bond) return;

	_strategy->addParameter(&*bond, Bond::getMagicPhi,
	                        Bond::setMagicPhi, range, interval,
	"ph" + bond->shortDesc());
	_strategy->addParameter(&*bond, Bond::getMagicPsi,
	                        Bond::setMagicPsi, range, interval,
	"ps" + bond->shortDesc());

}

void Sampler::addBendAngle(BondPtr bond, double range, double interval)
{
	//    double number = fabs(range / interval);

	if (!bond->getRefineBondAngle())
	{
		return;
	}

	_strategy->addParameter(&*bond, Bond::getBendAngle, Bond::setBendAngle,
	                        range, interval, "b0" + bond->shortDesc());
	_strategy->addParameter(&*bond, Bond::getCirclePortion, Bond::setCirclePortion,
	                        range, interval, "b1" + bond->shortDesc());

	_bonds.push_back(bond);
}

void Sampler::addAbsolutePosition(AbsolutePtr abs, double range, double interval)
{
	if (abs->getClassName() != "Absolute")
	{
		return;
	}

	_strategy->addParameter(&*abs, Absolute::getPosX, Absolute::setPosX,
	                        range, interval, "pos_x");
	_strategy->addParameter(&*abs, Absolute::getPosY, Absolute::setPosY,
	                        range, interval, "pos_y");
	_strategy->addParameter(&*abs, Absolute::getPosZ, Absolute::setPosZ,
	                        range, interval, "pos_z");
}

void Sampler::addRotamer(Sidechain *side, double range, double interval)
{
	//    double number = fabs(range / interval);
	_strategy->addParameter(side, Sidechain::getRotamerExponent,
	                        Sidechain::setRotamerExponent,
	range, interval, "rot_exp");
}


void Sampler::addAbsoluteBFactor(AbsolutePtr abs, double range, double interval)
{
	if (abs->getClassName() != "Absolute")
	{
		return;
	}

	_strategy->addParameter(&*abs, Absolute::getB, Absolute::setB,
	                        range, interval, "bfactor");
}

void Sampler::addSampled(std::vector<AtomPtr> atoms)
{
	for (int i = 0; i < atoms.size(); i++)
	{
		addSampled(atoms[i]);
	}
}

void Sampler::addSampledAtoms(AtomGroupPtr group, std::string conformer)
{
	if (!group)
	{
		return;
	}

	for (int i = 0; i < group->atomCount(); i++)
	{
		if (conformer.length() &&
		    group->atom(i)->getAlternativeConformer() != conformer)
		{
			continue;
		}

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
		/* No repeats! */
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
	_crystal = crystal;
	_real2Frac = crystal->getReal2Frac();
	_fft = crystal->getFFT();
}

void Sampler::setupCloseAtoms()
{
	_crystal->clearCloseCache();

	for (int i = 0; i < sampleSize(); i++)
	{
		_crystal->getCloseAtoms(_sampled[i], 9.0, true);
	}
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
		int paramCount = _strategy->parameterCount();
		int cycles = paramCount / 2;
		if (cycles < 10) cycles = 10;

		_strategy->setCycles(cycles);
	}

	int paramCount = _strategy->parameterCount();

	if (!_silent)
	{
		std::cout << "Sampling " << sampleSize() << " atoms, ";
		std::cout << "refining " << paramCount << " parameters." << std::endl;
	}

	if (_scoreType == ScoreTypeCorrel)
	{
		int cycles = 16 + paramCount;

		_strategy->setCycles(cycles);

		setupCloseAtoms();
		AtomGroup::scoreWithMapGeneral(_scoreType, _crystal, true, _sampled);
	}

	if (sampleSize() && _strategy->parameterCount())
	{
		_strategy->setJobName(_jobName);
		_strategy->refine();
	}

	_silent = false;
	_scoreType = ScoreTypeCorrel;

	bool changed = _strategy->didChange();

	if (clear)
	{
		_strategy = RefinementStrategyPtr();
	}

	_bonds.clear();
	_sampled.clear();
	_unsampled.clear();

	return changed;
}

double Sampler::getScore()
{
	if (!_sampled.size())
	{
		return 0;
	}

	if (_scoreType == ScoreTypeModelPos || _scoreType == ScoreTypeModelRMSDZero)
	{
		double score = 0;
		double count = 0;

		for (int i = 0; i < sampleSize(); i++)
		{
			double oneScore = 0;
			
			switch (_scoreType)
			{
				case ScoreTypeModelPos:
				oneScore = _sampled[i]->posDisplacement();
				break;
				
				case ScoreTypeModelRMSDZero:
				oneScore = _sampled[i]->fullPositionDisplacement();
				break;
				
				default:
				break;
			}

			score += oneScore;
			count += 1;
		}

		return score / count;
	}

	return AtomGroup::scoreWithMapGeneral(_scoreType, _crystal, false, _sampled);
}


