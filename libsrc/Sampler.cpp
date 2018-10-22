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
#include "RefinementStepSearch.h"
#include "Bond.h"
#include "Atom.h"
#include "Crystal.h"
#include "Element.h"
#include "FileReader.h"
#include "Absolute.h"
#include "Polymer.h"
#include "Monomer.h"
#include "Backbone.h"

/** \endcond */

Sampler::Sampler()
{
	_mock = false;
	_scoreType = ScoreTypeCorrel;
	_refinedMagicAxisCount = 0;
	_silent = false;
}

void Sampler::addAtomsForBond(BondPtr bond)
{
	addSampled(bond->getMinor());
	addSampled(bond->getMajor());

	int count = 0;

	/* Stop if the bond has no parent */
	if (!bond->getParentModel()->isBond())
	{
		return;
	}

	BondPtr parent = ToBondPtr(bond->getParentModel());
	int group = -1;
	int num = parent->downstreamBondNum(&*bond, &group);

	/* if it's the oldest sibling, we need to add all the
	 * siblings as this torsion controls the others */
	if (num == 0)
	{
		for (int j = 1; j < parent->downstreamBondCount(group); j++)
		{
			AtomPtr downAtom = parent->downstreamAtom(group, j);
			
			addSampled(downAtom);
			count++;
		}
	}

	for (int j = 0; j < bond->extraTorsionSampleCount(); j++)
	{
		addSampled(bond->extraTorsionSample(j));
		count++;
	}
	
	if (bond->downstreamBondGroupCount())
	{
		for (size_t j = 0; j < bond->downstreamBondGroupCount(); j++)
		{
			for (size_t i = 0; i < bond->downstreamBondCount(j); i++)
			{
				AtomPtr downAtom = bond->downstreamAtom(j, i);
				addSampled(downAtom);
				count++;
			}
		}
	}
}

void Sampler::addCustomParameter(void *object, Getter getter, Setter setter,
                                 double range, double interval,
                                 std::string name)
{
	_strategy->addParameter(object, getter, setter, range, interval, name);
}

int Sampler::hasParameter(ParamOptionType type)
{
	return (_params.count(type));
}

void Sampler::addParamsForBond(BondPtr bond, bool even)
{
	int mult = (even ? 1 : -1);
	for (ParamMap::iterator it = _params.begin(); it != _params.end(); it++)
	{
		ParamOptionType option = it->first;
		double range = it->second;
		
		if (range < 1e-6)
		{
			continue;
		}

		switch (option)
		{
			case ParamOptionTorsion:
			addTorsion(bond, deg2rad(range) * mult, deg2rad(0.005));
			break;

			case ParamOptionBondAngle:
			addBendAngle(bond, deg2rad(range) * mult, deg2rad(0.005));
			break;

			case ParamOptionKick:
			addKick(bond, range * mult, 0.001);
			break;

			case ParamOptionMagicAngles:
			addMagicAngle(bond, deg2rad(range) * mult, deg2rad(1.0));
			break;

			default:
			break;
		}
	}
}

BondPtr Sampler::setupThoroughSet(BondPtr bond, bool addBranches)
{
	reportInDegrees();
	setScoreType(ScoreTypeCorrel);

	setJobName("torsion_set_" + bond->shortDesc());

	int bondCount = 0;
	BondPtr topBond = BondPtr();

	if (_params.count(ParamOptionNumBonds) == 0)
	{
		return BondPtr();
	}
	
	int bondNum = _params[ParamOptionNumBonds] + 0.2;

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
		addAtomsForBond(bond);
		
		if (!bond->isRefinable())
		{
			/* No hope! Give up! */
			continue;
		}

		addParamsForBond(bond, (num % 2));

		for (int j = 0; j < bond->downstreamBondGroupCount(); j++)
		{
			/* Take the chosen group and check for futures */
			for (int i = 0; i < bond->downstreamBondCount(j); i++)
			{
				AtomPtr downstreamAtom = bond->downstreamAtom(j, i);
				BondPtr nextBond = ToBondPtr(downstreamAtom->getModel());
				
				if (!topBond)
				{
					topBond = nextBond;
				}

				BondInt entry;
				entry.bond = nextBond;
				entry.num = num - 1;
				
				if (addBranches || (!addBranches && i == 0))
				{
					remaining.push_back(entry);
				}
			}
		}
	}
	
	if (topBond && !topBond->isRefinable())
	{
		return BondPtr();
	}
	
	return topBond;
}

void Sampler::setupGrid()
{
	_strategy = RefinementStrategyPtr(new RefinementGridSearch());
	_strategy->setEvaluationFunction(Sampler::score, this);
}

void Sampler::setupNelderMead()
{
	_strategy = RefinementStrategyPtr(new RefinementNelderMead());
	_strategy->setEvaluationFunction(Sampler::score, this);
	_strategy->setCycles(20);
}

void Sampler::setupStepSearch()
{
	_strategy = RefinementStrategyPtr(new RefinementStepSearch());
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

	_strategy->addParameter(&*bond, Bond::getTorsion, 
	                        Bond::setTorsion,
	                        range, interval,
	"t" + bond->shortDesc());

	_bonds.push_back(bond);
}

void Sampler::addKick(BondPtr bond, double range, double interval)
{
	if (!bond) return;
	
	if (!bond->getRefineFlexibility()) return;

	_strategy->addParameter(&*bond, Bond::getKick, Bond::setKick,
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

void Sampler::addMagicAngle(BondPtr bond, double range, double interval)
{
	if (!bond) return;
	if (!bond->getRefineFlexibility()) return;
	
	if (Bond::getKick(&*bond) < 1e-6 && 
	    (!hasParameter(ParamOptionKick) || _params[ParamOptionKick] < 1e-6))
	{
		return;
	}
	
	_strategy->addParameter(&*bond, Bond::getMagicPhi,
	                        Bond::setMagicPhi, range, interval,
	"ph" + bond->shortDesc());
	_strategy->addParameter(&*bond, Bond::getMagicPsi,
	                        Bond::setMagicPsi, range, interval,
	"ps" + bond->shortDesc());

}

void Sampler::addBendAngle(BondPtr bond, double range, double interval)
{
	if (!bond->getRefineBondAngle() || bond->isFixed())
	{
		return;
	}

	_strategy->addParameter(&*bond, Bond::getBendAngle, Bond::setBendAngle,
	                        range, interval, "b0_" + bond->shortDesc());
	
	ModelPtr parent = bond->getParentModel();

	if (parent && parent->isBond())
	{
		BondPtr pBond = ToBondPtr(parent);
		int result = pBond->downstreamBondNum(&*bond, NULL);

		if (result > 0)
		{
			_strategy->addParameter(&*bond, Bond::getCirclePortion, 
			                        Bond::setCirclePortion,
			range, interval, "b1_" + bond->shortDesc());
		}
	}

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

void Sampler::setupScoreWithMap()
{
	_workspace.scoreType = _scoreType;
	_workspace.crystal = _crystal;
	_workspace.selectAtoms = _sampled;
	_workspace.segment = FFTPtr();
	_workspace.ave = empty_vec3();
	_workspace.basis = make_mat3x3();
	_workspace.flag = MapScoreFlagNone;
		
	AtomGroup::scoreWithMapGeneral(&_workspace);
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
		std::cout << "Refining " << paramCount << " parameters." << std::endl;
		std::cout << "Sampling " << sampleSize() << " atoms: ";
		
		for (int i = 0; i < _sampled.size(); i++)
		{
			std::cout << _sampled[i]->shortDesc() << ", ";
		}
		
		std::cout << std::endl;
	}

	if (_scoreType == ScoreTypeCorrel ||
	    _scoreType == ScoreTypeRFactor ||
	    _scoreType == ScoreTypeMultiply ||
		_scoreType == ScoreTypeHappiness)
	{
		int cycles = 16 + paramCount;

		_strategy->setCycles(cycles);
		setupCloseAtoms();
		setupScoreWithMap();
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
	_crystal->clearCloseCache();

	return changed;
}

double Sampler::getScore()
{
	if (!_sampled.size() || _scoreType == ScoreTypeZero)
	{
		return 0;
	}
	
	if (_scoreType == ScoreTypeCentroid)
	{
		vec3 orig_cent = empty_vec3();
		vec3 new_cent = empty_vec3();
		double count = 0;
		
		for (int i = 0; i < _bonds.size() && i < 1; i++)
		{
			_bonds[i]->propagateChange(-1, true);
		}

		for (int i = 0; i < sampleSize(); i++)
		{
			vec3 init = _sampled[i]->getPDBPosition();
			_sampled[i]->getModel()->refreshPositions();
			vec3 pos = _sampled[i]->getAbsolutePosition();
			
			vec3_add_to_vec3(&orig_cent, init);
			vec3_add_to_vec3(&new_cent, pos);
			count++;
		}

		vec3_mult(&orig_cent, 1 / count);
		vec3_mult(&new_cent, 1 / count);
		
		vec3 diff = vec3_subtract_vec3(orig_cent, new_cent);

		return vec3_sqlength(diff);
	}
	
	if (_scoreType != ScoreTypeCorrel 
	    && _scoreType != ScoreTypeRFactor
		&& _scoreType != ScoreTypeHappiness)
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
				oneScore = _sampled[i]->getBFactor();
				break;
				
				case ScoreTypeBFactorAgreement:
				oneScore = pow(_sampled[i]->getBFactor() -
				            _sampled[i]->getInitialBFactor() + 6.3, 2.0);
				break;
				
				case ScoreTypeRMSDZero:
				oneScore = _sampled[i]->getModel()->smallness();
				break;
				
				default:
				break;
			}

			score += oneScore;
			count += 1;
		}
		
		return score / count;
	}

	double score = AtomGroup::scoreWithMapGeneral(&_workspace);
	
	return score;
}


