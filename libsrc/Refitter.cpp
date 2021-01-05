//
//  Reflex.cpp
//  vagabond
//
//  Created by Helen Ginn, 2018
//  Copyright © 2018 Helen Ginn. All rights reserved.
//

#include "Refitter.h"
#include <hcsrc/RefinementGridSearch.h>
#include "Bond.h"
#include "Monomer.h"
#include "Polymer.h"
#include "Backbone.h"
#include "Atom.h"
#include "Options.h"
#include "ParamBand.h"
#include "Shouter.h"
#include <hcsrc/Timer.h>

double ramachandran_angles[] = {
	-180, 180,
	-135, 180,
	-90, 180,
	-45, 180,
	-180, 135,
	-135, 135,
	-90, 135,
	-45, 135,
	-180, 90,
	-135, 90,
	-90, 90,
	-45, 90,
	-135, 45,
	-90, 45,
	-45, 45,
	-135, 0,
	-90, 0,
	-45, 0,
	-135, 45,
	-90, 45,
	-45, 45,
	0, 45,
	-135, -90,
	-90, -90,
	-45, -90,
	0, -90,
	-135, -135,
	-90, -135,
	-45, -135,
	45, 90,
	45, 45,
	45, 0,
	55, -160,
};

Refitter::Refitter(BondPtr bond, BondPtr target, bool forwards)
{
	_bond = bond;
	_mStart = _bond->getMinor()->getResidueNum();
	_group = _bond->makeAtomGroup();
	
	_forwards = forwards;
	_mEnd = _mStart;

	for (int i = 0; i < _group->atomCount(); i++)
	{
		AtomPtr atom = _group->atom(i);
		int mon = atom->getResidueNum();
		
		if ((forwards && mon > _mEnd) || (!forwards && mon < _mEnd))
		{
			_mEnd = mon;
		}
	}
	
	_mLand = (_forwards ? _mEnd - 1 : _mEnd + 1);
	_polymer = _bond->getAtom()->getMonomer()->getPolymer();
	
	if (fabs(_mStart - _mEnd) < 3)
	{
		std::cout << "Only fitting less than three residues?" << std::endl;
	}
	
	std::cout << "Residues " << _mStart << " to " << _mEnd <<
	" and a landing residue of " << _mLand << std::endl;
	
	_target = target;
	_degstep = 45;
	_pepAtom = (_forwards ? "N" : "C");

	vec3 pos = _group->centroid();
//	Options::getRuntimeOptions()->focusOnPosition(pos);
}

double Refitter::accumulativeBondLength(BondPtr current)
{
	double sum = 0;
	
	/* Return whatever is not the next Ramachandran angle */
	while (current)
	{
		if (!current->downstreamBondGroupCount() ||
		    !current->downstreamBondCount(0))
		{
			break;
		}
		
		sum += Bond::getBondLength(&*current);
		current = current->downstreamBond(0, 0);
	}
	
	return sum;
}

void Refitter::strictPrune(double angle)
{
	size_t purged = 0;
	AtomPtr last = _group->findAtoms("CA", _mEnd)[0];
	AtomPtr prev = _group->findAtoms((_forwards ? "O" : "N"), _mEnd)[0];
	vec3 init = last->getInitialPosition();
	vec3 previnit = prev->getInitialPosition();
	vec3 diff = vec3_subtract_vec3(previnit, init);

	for (int i = 0; i < _trialList.size(); i++)
	{
		ParamBandPtr trial = _trialList[i].list;
		applyTrial(trial);

		last->getModel()->refreshPositions();
		prev->getModel()->refreshPositions();
		vec3 pos = last->getAbsolutePosition();
		vec3 prevpos = prev->getAbsolutePosition();
		
		vec3 mydiff = vec3_subtract_vec3(prevpos, pos);
		double diffang = vec3_angle_with_vec3(diff, mydiff);

		if (diffang > deg2rad(angle))
		{
			purged++;
			_trialList.erase(_trialList.begin() + i);
			i--;
		}
	}
	
	std::cout << "Strictly purged " << purged << ", leaving "
	<< _trialList.size() << std::endl;

}

bool trial_lt_trial(const BondTrial a, const BondTrial b)
{
	return (a.score < b.score);
}

void Refitter::pruneUnreasonable(double distance, bool landing)
{
	int start = landing ? _mEnd : _mStart;
	int ending = _mEnd;
	int dir = _forwards ? 1 : -1;
	ending += dir;
	
	std::cout << "Pruning unreasonable trials between " << start
	<< " and " << ending << ", direction " << dir << std::endl;

	for (int i = 0; i < _trialList.size(); i++)
	{
		ParamBandPtr trial = _trialList[i].list;
		applyTrial(trial);
		
		double count = 0;
		double sum = 0;
		
		for (int j = start; j != ending; j += dir)
		{
			AtomList checks = _group->findAtoms("CA", j);

			if (!checks.size())
			{
				std::cout << "Warning: _group gave no atoms " << j << std::endl;
			}

			AtomPtr check = checks[0];
			vec3 aim = check->getInitialPosition();
			vec3 pos = check->getAbsolutePosition();
			vec3 diff = vec3_subtract_vec3(pos, aim);

			double length = vec3_length(diff);
			sum += length;
			count++;
		}
		
		sum /= count;
		_trialList[i].score = sum;
	}

	std::sort(_trialList.begin(), _trialList.end(), trial_lt_trial);
	
	int remove = 80 * fabs(_mStart - _mLand);
	int removed = _trialList.size() - remove;
	if (removed < 0) removed = 0;

	if (distance > 0)
	{
		for (int i = 0; i < _trialList.size(); i++)
		{
			if (_trialList[i].score > distance)
			{
				_trialList.erase(_trialList.begin() + i);
				i--;
				removed++;
			}
		}
	}
	else if (removed > 0)
	{
		_trialList.erase(_trialList.begin() + remove, _trialList.end());
	}
	
	double furthest = _trialList.back().score;

	std::cout << "Removed " << removed << " over " << furthest <<  "." <<
	std::endl;
}

void Refitter::refit()
{
	makePeptides180();
	
	BondPtr current = _bond;

	if (_forwards)
	{
		current = findNextBond(current);
	}
	
	int res = _mStart;
	
	/* First one needs to be started with no seeds */
	std::vector<BondTrial> newTrials;
	current = generateCandidates(ParamBandPtr(), res, &newTrials);
	res += _forwards ? 1 : -1;
	std::cout << "First trials: " << newTrials.size() << std::endl;
	_trialList = newTrials;
	
	while (res != _mLand)
	{
		std::cout << "Moving to residue " << res << std::endl;
		std::vector<BondTrial> allTrials;
		BondPtr next = current;

		for (int i = 0; i < _trialList.size(); i++)
		{
			next = generateCandidates(_trialList[i].list, res,
			                          &allTrials);
		}

		current = next;
		_trialList = allTrials;

		pruneUnreasonable();
		res += _forwards ? 1 : -1;

		std::cout << "Number of trials: " << _trialList.size() << std::endl;
		
		if (!_trialList.size())
		{
			break;	
		}
	}
	
	std::cout << "We have " << _trialList.size() << " surviving"
	" candidates." << std::endl;

	if (!_trialList.size())
	{
		return;	
	}

	Timer timer("centroid refinement", true);

	for (int i = 0; i < _trialList.size(); i++)
	{
		refineCentroid(&_trialList[i], 60);
	}
	
	timer.report();
	
//	removeDuplicates();
	pruneUnreasonable(1.5, true);
	strictPrune(40);
	
	std::cout << "We have " 
	<< _trialList.size() << " surviving candidates." << std::endl;
	
//	removeLanding();
	
	addSolutions();

	BondPtr parent = ToBondPtr(_bond->getParentModel());
	parent->equaliseOccupancies();
}

/*
void Refitter::removeEverything()
{
	Options::pauseGUIFishing(true);

	AtomList list = _polymer->findAtoms("C", _mEnd);
	if (list.size())
	{
		ToBondPtr(list[0]->getModel())->setResetOccupancy(false);
	}
	
	BondPtr parent = ToBondPtr(_bond->getParentModel());
	_bond->destroy();
	AtomList orig = _polymer->findAtoms(pAtom->getAtomName(),
	                                    pAtom->getResidueNum());
	
	if (orig.size())
	{
		ToBondPtr(orig[0]->getModel())->setResetOccupancy(true);
	}

	Options::pauseGUIFishing(false);

}
*/

void Refitter::removeLanding()
{
	Options::pauseGUIFishing(true);

	AtomList list = _polymer->findAtoms(_pepAtom, _mEnd);
	
	/* Remove the resetting of occupancy for primary conformer */
	if (list.size())
	{
		ToBondPtr(list[0]->getModel())->setResetOccupancy(false);
	}
	
	AtomList landings = _group->findAtoms("C", _mLand);

	AtomPtr land = landings[0];
	BondPtr bond = ToBondPtr(land->getModel());
	AtomPtr pAtom = bond->getMajor();
	bond->destroy();
	AtomList orig = _polymer->findAtoms(pAtom->getAtomName(),
	                                    pAtom->getResidueNum());
	
	if (orig.size())
	{
		ToBondPtr(orig[0]->getModel())->setResetOccupancy(true);
	}

	Options::pauseGUIFishing(false);
}

void Refitter::addSolutions()
{
	if (_trialList.size() == 0)
	{
		return;
	}

	applyTrial(_trialList[0].list);
	
	for (int i = 1; i < _trialList.size(); i++)
	{
		_bond->splitBond();
		applyTrial(_trialList[i].list, false, i);
	}
	
	_bond->getParentModel()->propagateChange(-1, true);
}

void Refitter::removeDuplicates()
{
	int count = 0;

	for (int i = 0; i < _trialList.size() - 1; i++)
	{
		ParamBandPtr trial = _trialList[i].list;
		for (int j = i + 1; j < _trialList.size(); j++)
		{
			ParamBandPtr other = _trialList[j].list;
			bool different = false;
			
			if (trial->objectCount() != other->objectCount())
			{
				std::cout << "!!!" << std::endl;
			}

			for (int k = 0; k < trial->objectCount(); k++)
			{
				double iTor = trial->baseValueForObject(k);
				double jTor = other->baseValueForObject(k);

				if (fabs(iTor - jTor) > deg2rad(_degstep) * 1.1)
				{
					different = true;
					break;
				}
			}

			if (!different)
			{
				/* one of these needs to die */
				bool i_bt_j = (_trialList[i].score < _trialList[j].score);
				int kill = (i_bt_j ? j : i);
				_trialList.erase(_trialList.begin() + kill);
				j--;
				
				if (!i_bt_j)
				{
					i--;
				}
				
				count++;
				break;
			}
		}
	}

	std::cout << "Removed " << count << " duplicate solutions." << std::endl;
}

void Refitter::refineCentroid(BondTrial *trial, int tries)
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	applyTrial(trial->list); 
		
	/* stage 0 - refine centroid of backbone
	 * stage 1 - refine centroid of sidechain
	 * stage 2 - refine per bond */
	int stage = 0;

	BondPtr last = trial->last;
	AtomPtr lastAtom = last->getMinor();
	MonomerPtr monomer = lastAtom->getMonomer();
	lastAtom = _group->findAtoms("CA", _mEnd)[0];
	monomer = lastAtom->getMonomer();
	BackbonePtr bone = monomer->getBackbone();
	SidechainPtr side = monomer->getSidechain();
	AtomGroupPtr add = bone; 
	
	for (int i = 0; i < bone->atomCount(); i++)
	{
		bone->atom(i)->setFromPDB(true);
	}
	
	double multiplier = 8;
	
	for (int i = 0; i < tries && (stage < 2); i++)
	{
		setCrystal(crystal);

		ScoreType type = (stage > 1) ? ScoreTypeModelPos : ScoreTypeCentroid;
		
		setScoreType(type);
		setupNelderMead();
		setCycles(100);
		reportInDegrees();
		setSilent();
		addSampledAtoms(add, "b");
		
		double angle = getScore() * multiplier;
		double limit = angle / 100;
		if (angle < 1) angle = 1;
		if (limit < 0.1) limit = 0.1;
		angle *= rand() % 2 ? 1 : -1;

		int start = _mStart;
		int end = _mEnd;
		
		for (int j = start; j != end; j += (_forwards ? 1 : -1))
		{
			AtomList nitro = _group->findAtoms("N", j);
			AtomList carbon = _group->findAtoms("C", j);

			BondPtr mf = ToBondPtr(nitro[0]->getModel());
			BondPtr ml = ToBondPtr(carbon[0]->getModel());
			
			addTorsion(mf, deg2rad(angle), deg2rad(limit));
			addTorsion(ml, deg2rad(angle), deg2rad(limit));
		}

		for (int j = 0; j < bone->atomCount() && (stage > 1); j++)
		{
			AtomPtr atom = bone->atom(j);
			BondPtr bond = ToBondPtr(atom->getModel());
			
			addTorsion(bond, deg2rad(angle), deg2rad(limit));
		}
		
		bool changed = sample();
		
		if (!changed && (stage == 1))
		{
			/* we need to go up to per-bond */
			stage = 2;
		}
		else if ((!changed || (i == tries - 1)) && (stage == 0))
		{
			/* we need to switch to side chain */
			stage = 1;
			multiplier = 4;
			i = 0;
			add = AtomGroupPtr(new AtomGroup());
			
			AtomList atoms = _group->findAtoms("CB", 
			                                   monomer->getResidueNum());
			
			if (atoms.size())
			{
				add->addAtom(atoms[0]);
			}

			atoms = _group->findAtoms("CA", monomer->getResidueNum());
			if (atoms.size())
			{
				add->addAtom(atoms[0]);
			}
		}
	}

	collectTorsionsIntoTrial(trial);

	AtomList landings = _group->findAtoms("CA", _mLand);
	
	if (!landings.size())
	{
		shout_at_helen("No landing atom for some reason.");
	}

	AtomPtr landing = landings[0];
	BondPtr landBond = ToBondPtr(landing->getModel());
	AtomGroupPtr all = landBond->makeAtomGroup(last);

	for (int i = 0; i < all->atomCount(); i++)
	{
		all->atom(i)->setFromPDB(true);
	}
	
	all->addParamType(ParamOptionTorsion, 2.0);
	all->addParamType(ParamOptionNumBonds, 5);
	all->addParamType(ParamOptionMaxTries, tries);

	all->refine(crystal, RefinementModelPos);
	collectTorsionsIntoTrial(trial);
	
	_solutions.push_back(collectTorsionsIntoSolution());
	
	for (int i = 0; i < all->atomCount(); i++)
	{
		all->atom(i)->setFromPDB(false);
	}
}

ParamBandPtr Refitter::collectTorsionsIntoSolution()
{
	ParamBandPtr band = ParamBandPtr(new ParamBand());
	band->setPrivateGetter(Bond::getTorsion);
	band->setPrivateSetter(Bond::setTorsion);
	
	for (int i = 0; i < _group->atomCount(); i++)
	{
		band->addObject(&*_group->atom(i)->getModel(), 1.);
	}
	
	band->prepare();
	return band;
}

void Refitter::collectTorsionsIntoTrial(BondTrial *trial)
{
	trial->list->prepare();
}

BondPtr Refitter::findNextBond(BondPtr bond)
{
	BondPtr current = bond;
	
	/* Return whatever is not the next Ramachandran angle */
	while (current)
	{
		if (current != bond && 
		    !(current->connectsAtom("CA") && current->connectsAtom(_pepAtom)))
		{
			break;
		}
		
		if (!current->downstreamBondGroupCount() ||
		    !current->downstreamBondCount(0))
		{
			return BondPtr();
		}
		
		current = current->downstreamBond(0, 0);
	}
	
	return current;
}

typedef struct
{
	size_t index;
	size_t count;	
	double value;
} IntPair;

inline bool ip_gt_ip(const IntPair &ip1, const IntPair &ip2)
{
	if (ip1.count == 0 && ip2.count > 0)
	{
		return false;
	}

	return ip1.value < ip2.value;
}

BondPtr Refitter::advanceBondConformer(BondPtr bond, int advance)
{
	AtomPtr minor = bond->getMinor();
	
	MonomerPtr monomer = minor->getMonomer();
	AtomList atoms = monomer->findAtoms(minor->getAtomName());
	
	if (advance + 1 >= atoms.size())
	{
		return BondPtr();
	}

	BondPtr next = ToBondPtr(atoms[advance + 1]->getModel());
	
	return next;
}

void Refitter::applyTrial(ParamBandPtr trial, bool propagate, int advance)
{
	int count = 0;

	if (trial)
	{
		for (int j = 0; j < trial->objectCount(); j++)
		{
			void *bObj = trial->object(j);
			BondPtr curr = static_cast<Bond *>(bObj)->shared_from_this();
			double val = trial->baseValueForObject(j);

			if (advance > 0)
			{
				curr = advanceBondConformer(curr, advance);
			}
			
			if (!curr)
			{
				continue;
			}

			Bond::setTorsion(&*curr, val);
			count++;
		}
	}

	if (propagate)
	{
		_bond->propagateChange(-1, true);
	}
}

BondPtr Refitter::generateCandidates(ParamBandPtr trial, int residue,
                           std::vector<BondTrial> *populate)
{
	AtomList nitro = _group->findAtoms("N", residue);
	AtomList carbon = _group->findAtoms("C", residue);
	
	if (!nitro.size() || !carbon.size())
	{
		return BondPtr();
	}
	
	ModelPtr mf = nitro[0]->getModel();
	ModelPtr ml = carbon[0]->getModel();
	BondPtr former = ToBondPtr(mf);
	BondPtr latter = ToBondPtr(ml);
	
	if (!former || !latter)
	{
		return BondPtr();
	}

	double fTmp = Bond::getTorsion(&*former);
	double lTmp = Bond::getTorsion(&*latter);
	
	if (trial)
	{
		applyTrial(trial, false);
	}
	
	size_t size = sizeof(ramachandran_angles) / sizeof(double);
	
	for (size_t i = 0; i < size; i += 2)
	{
		ParamBandPtr newTrial = ParamBandPtr(new ParamBand());
		newTrial->setPrivateGetter(Bond::getTorsion);
		newTrial->setPrivateSetter(Bond::setTorsion);

		if (trial)
		{
			newTrial = ParamBandPtr(new ParamBand(*trial));
		}

		newTrial->addObject(&*former, 1);
		newTrial->addObject(&*latter, 1);

		double iDeg = ramachandran_angles[i];
		double jDeg = ramachandran_angles[i + 1];

		Bond::setTorsion(&*latter, deg2rad(iDeg));
		Bond::setTorsion(&*former, deg2rad(jDeg));

		newTrial->prepare();

		if (populate)
		{
			BondTrial trial;
			trial.score = 0;
			trial.list = newTrial;
			trial.last = _forwards ? latter : former; 
			populate->push_back(trial);
		}
	}

	return BondPtr();
}

void Refitter::makePeptides180()
{
	std::cout << "Making all peptides 180 degrees." << std::endl;
	BondPtr current = _bond;
	
	while (current)
	{
		if (current->connectsAtom("CA") && current->connectsAtom(_pepAtom))
		{
			Bond::setTorsion(&*current, deg2rad(180));
			current->setFixed(true);
		}
		
		if (!current->downstreamBondGroupCount() ||
		    !current->downstreamBondCount(0))
		{
			break;
		}
		
		current = current->downstreamBond(0, 0);
	}
	
	_bond->propagateChange(-1, true);
	
	std::cout << "Done." << std::endl;
}
