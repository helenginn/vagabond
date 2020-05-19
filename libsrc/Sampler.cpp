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
#include "SVDBond.h"
#include "Atom.h"
#include "Crystal.h"
#include "Element.h"
#include "FileReader.h"
#include "Absolute.h"
#include "Anchor.h"
#include "Polymer.h"
#include "Monomer.h"
#include "Backbone.h"
#include "Twist.h"
#include "Balance.h"

#include <iomanip>

/** \endcond */

Sampler::Sampler()
{
	_svd = NULL;
	_shouldSave = false;
	_begin = 0;
	_cycles = 0;
	_verbose = false;
	_improv = 0;
	_changed = false;
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
	
	/* Downstream bonds important if not backwards */
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
	
}

void Sampler::addTwistShift(ExplicitModelPtr eModel, AtomGroupPtr clearGroup)
{
	eModel->addShifts(_strategy, clearGroup);
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
		
		if (fabs(range) < 1e-6)
		{
			continue;
		}

		double degrange = deg2rad(range);
		double tol = degrange / 100.;
		
		if (tol < deg2rad(0.005))
		{
			tol = deg2rad(0.005);
		}

		switch (option)
		{
			case ParamOptionTorsion:
			addTorsion(bond, degrange * mult, tol);
			break;

			case ParamOptionTwist:
			addTwist(bond, degrange * mult, tol);
			break;

			case ParamOptionBondAngle:
			addBendAngle(bond, degrange * mult, tol);
			break;

			case ParamOptionCirclePortion:
			addBendAngle(bond, degrange * mult, tol, true);
			break;

			case ParamOptionKick:
			addKick(bond, range * mult, 0.001);
			break;

			case ParamOptionMagicAngles:
			addMagicAngle(bond, range * mult, deg2rad(1.0));
			break;

			default:
			break;
		}
	}
}

BondPtr Sampler::setupThoroughSet(BondPtr fbond, bool addBranches)
{
	reportInDegrees();

	setJobName("torsion_set_" + fbond->shortDesc());

	int bondCount = 0;
	BondPtr topBond = BondPtr();

	if (_params.count(ParamOptionNumBonds) == 0)
	{
		return BondPtr();
	}
	
	int bondNum = _params[ParamOptionNumBonds] + 0.2;
	int atomNum = bondNum;
	if (addBranches)
	{
		int extraNum = _params[ParamOptionExtraAtoms];
		atomNum = _params[ParamOptionNumBonds] + extraNum;
	}

	/* Logic: bond add map keeps remaining bonds (and remaining
 * 	additions) in memory for adding (allowing branches). When a
 * 	branched bond needs adding it is appended to this list with
 * 	one fewer remaining additions and the old entry deleted. */

	std::vector<BondInt> remaining;
	BondInt entry;
	entry.bond = fbond;
	entry.bondNum = bondNum;
	entry.atomNum = atomNum;
	remaining.push_back(entry);
	addSampled(fbond->getMajor());

	while (remaining.size())
	{
		BondInt first = remaining[0];
		BondPtr bond = first.bond;
		int aNum = first.atomNum;
		int bNum = first.atomNum;

		remaining.erase(remaining.begin());

		if (aNum <= 0)
		{
			continue;	
		}
		
		bondCount++;
		
		addAtomsForBond(bond);
		checkOccupancyAndAdd(bond);
		
		if (!bond->isRefinable() && bond != fbond)
		{
			/* No hope! Give up! Unless first bond */
			continue;
		}

		if (bond->isRefinable() && entry.bondNum >= 0)
		{
			addParamsForBond(bond, (bNum % 2));
		}

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
				entry.bondNum = bNum - 1;
				entry.atomNum = aNum - 1;
				
				/* exception: carbonyl oxygen we add branch anyway */
				bool isCarbonyl = bond->getMinor()->getAtomName() == "C";
				
				if (addBranches || (!addBranches && i == 0) || isCarbonyl)
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
	ToGridPtr(_strategy)->setWritePNG();
	ToGridPtr(_strategy)->setWriteCSV();
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

void Sampler::checkOccupancyAndAdd(BondPtr bond)
{
	BondPtr parent = bond;
	int count = 0;

	while (parent->getParentModel()->isBond()
	    && hasParameter(ParamOptionOccupancy) && count < 5)
	{
		parent = ToBondPtr(parent->getParentModel());

		if (parent->downstreamBondGroupCount() > 1)
		{
			bool clear = true;
			for (int i = 0; i < _balances.size(); i++)
			{
				if (_balances[i]->isFromBond(parent))
				{
					clear = false;
					break;
				}
			}

			if (clear)
			{
				BalancePtr balance = BalancePtr(new Balance(parent));
				balance->addParamsToStrategy(_strategy);
				_balances.push_back(balance);
				break;
			}
		}
		
		count++;
	}
}

void Sampler::addOccupancy(BondPtr bond, double range, double interval)
{
	//    double number = fabs(range / interval);
	_strategy->addParameter(&*bond, Bond::getOccupancy, Bond::setOccupancy,
	                        range, interval,
	"occupancy");

	_bonds.push_back(bond);
}

void Sampler::addTwist(BondPtr bond, double range, double interval)
{
	if (bond->hasTwist())
	{
		TwistPtr twist = bond->getTwist();
		_strategy->addParameter(&*twist, Twist::getTwist,
		                        Twist::setTwist, range, interval,
		                        "tw" + bond->shortDesc());
		
		if (hasParameter(ParamOptionShift))
		{
			ExplicitModelPtr eModel = twist->getAppliedModel();
			addTwistShift(eModel, AtomGroupPtr());
			ParamMap::iterator it = _params.find(ParamOptionShift);
			_params.erase(it);
		}
	}
	else
	{
		addTorsion(bond, range, interval);
	}
}

void Sampler::addTorsion(BondPtr bond, double range, double interval)
{
	if (!bond || !bond->isBond())
	{
		return;
	}

	if (!bond->isRefinable() | !bond->isUsingTorsion())
	{
		if (!_silent)
		{
			std::cout << "Can't refine torsion for " << bond->shortDesc() << std::endl;
		}
		return;
	}
	
	if (!_silent)
	{
		std::cout << "Adding torsion for " << bond->shortDesc() << std::endl;
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
	
	if (bond->hasWhack())
	{
		WhackPtr whack = bond->getWhack();
		_strategy->addParameter(&*whack, Whack::getKick, Whack::setKick,
		                        range, interval, "k" + bond->shortDesc());
	}
	else
	{
		_strategy->addParameter(&*bond, Bond::getKick, Bond::setKick,
		                        range, interval, "k" + bond->shortDesc());
	}


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

void Sampler::addBendAngle(BondPtr bond, double range, double interval,
                           bool circlePortionOnly)
{
	if (!bond->getRefineBondAngle())
	{
		return;
	}
	
	if (!circlePortionOnly)
	{
		_strategy->addParameter(&*bond, Bond::getBendAngle, Bond::setBendAngle,
		                        range, interval, "b0_" + bond->shortDesc());
	}
	
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

void Sampler::addAnchorParams(AnchorPtr anch)
{
	double range = 0.02;
	double interval = 0.002;
	_strategy->addParameter(&*anch, Anchor::getPosX, Anchor::setPosX,
	                        range, interval, "pos_x");
	_strategy->addParameter(&*anch, Anchor::getPosY, Anchor::setPosY,
	                        range, interval, "pos_y");
	_strategy->addParameter(&*anch, Anchor::getPosZ, Anchor::setPosZ,
	                        range, interval, "pos_z");
	
	range = deg2rad(1.0);
	interval = deg2rad(0.005);
	_strategy->addParameter(&*anch, Anchor::getAlpha, Anchor::setAlpha,
	                        range, interval, "alpha");
	_strategy->addParameter(&*anch, Anchor::getBeta, Anchor::setBeta,
	                        range, interval, "beta");
	_strategy->addParameter(&*anch, Anchor::getGamma, Anchor::setGamma,
	                        range, interval, "gamma");

	BondPtr nBond = ToBondPtr(anch->getNAtom()->getModel());
	BondPtr cBond = ToBondPtr(anch->getCAtom()->getModel());
	
	addAtomsForBond(nBond);
	addAtomsForBond(cBond);
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
	
	if (atom->getElectronCount() == 1 && _scoreType == ScoreTypeSavedPos
	|| atom->getElectronCount() == 1 && _scoreType == ScoreTypeModelPos)
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
//	_fft = crystal->getFFT();
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
	AtomGroupPtr sampled = AtomGroupPtr(new AtomGroup());
	
	for (int i = 0; i < sampleSize(); i++)
	{
		sampled->addAtom(_sampled[i]);
	}
	
	for (int i = 0; i < _includeForRefine.size(); i++)
	{
		sampled->addAtomsFrom(_includeForRefine[i]);
	}
	
	double sampling = getParameter(ParamOptionProteinSampling);
	sampled->addParamType(ParamOptionProteinSampling, sampling);
	copyParams(sampled);
	
	setup_space(&_workspace);
	_workspace.scoreType = _scoreType;
	_workspace.crystal = _crystal;
	_workspace.selectAtoms = sampled;
		
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
		if (cycles < 30) cycles = 30;
		_strategy->setCycles(cycles);
	}
	
	if (_params.count(ParamOptionCycles) > 0)
	{
		_strategy->setCycles(_params[ParamOptionCycles]);
	}

	if (_params.count(ParamOptionSVD) > 0)
	{
		_strategy->clearParameters();
		_svd = new SVDBond(_bonds, _sampled);
		_svd->setDoTorsion(true);
		_svd->setSilent(_silent);
		_svd->performSVD();
		double t = 0.2;
		if (_params.count(ParamOptionTorsion))
		{
			t = _params[ParamOptionTorsion];
		}
		_svd->addToStrategy(_strategy, t, false);
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

		if (_cycles > 0)
		{
			cycles = _cycles;
		}

		if (_params.count(ParamOptionSVD) > 0)
		{
			cycles *= 5;
		}

		_strategy->setCycles(cycles);
//		setupCloseAtoms();
		setupScoreWithMap();
	}
	
	if (_verbose)
	{
		_strategy->setVerbose(true);
	}
	
	if (_shouldSave)
	{
		_begin = getScore();
		_shouldSave = false;
	}

	if (sampleSize() && _strategy->parameterCount())
	{
		_strategy->setJobName(_jobName);

		int maxTries = getParameter(ParamOptionMaxTries);
		for (int i = 0; i < maxTries; i++)
		{
			_strategy->refine();
			if (!_strategy->didChange())
			{
				break;
			}
		}
	}

	double end = getScore();
	_improv = end - _begin;
	
	_silent = false;
	_scoreType = ScoreTypeCorrel;

	_changed = _strategy->didChange();
	

	if (clear)
	{
		_balances.clear();
		_strategy = RefinementStrategyPtr();
	}

	_cycles = 0;
	_bonds.clear();
	_sampled.clear();
	_crystal->clearCloseCache();
	delete _svd;
	_svd = NULL;

	return _changed;
}

void Sampler::saveScore()
{
	_shouldSave = true;
}

double Sampler::getScore()
{
	if (!_sampled.size() || _scoreType == ScoreTypeZero)
	{
		return 0;
	}
	
	if (_svd != NULL)
	{
		_svd->applyParameters();
	}

	for (int i = 0; i < _sampled.size(); i++)
	{
		AtomPtr s = _sampled[i];
		s->getModel()->propagateChange(0);
		s->getModel()->refreshPositions();
	}

	
	for (int i = 0; i < _balances.size(); i++)
	{
		_balances[i]->adjustment();
	}
	
	if (_scoreType == ScoreTypeCentroid)
	{
		vec3 orig_cent = empty_vec3();
		vec3 new_cent = empty_vec3();
		double count = 0;
		
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
				
				case ScoreTypeSavedPos:
				oneScore = _sampled[i]->posDisplacement(true);
				break;
				
				case ScoreTypeMouse:
				oneScore = _sampled[i]->posToMouse();
				break;
				
				default:
				break;
			}

			score += oneScore;
			count += 1;
		}
		
		return score / count;
	}

	double score = 0;
	if (sampleSize())
	{
		score = AtomGroup::scoreWithMapGeneral(&_workspace);
		score += constraint();
	}

	return score;
}

double Sampler::constraint()
{
	for (int i = 0; i < _bonds.size(); i++)
	{
		bool allowed = _bonds[i]->allowBondAngle();
		if (!allowed)
		{
			return FLT_MAX;
		}
	}

	return 0;
}
