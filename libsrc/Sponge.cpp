// Vagabond
// Copyright (C) 2019 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#include <algorithm>
#include "Sponge.h"
#include "Atom.h"
#include <hcsrc/Any.h>
#include "Options.h"
#include "ExplicitModel.h"
#include "WaterNetwork.h"

void Sponge::findConnections()
{
	CrystalPtr crystal = Options::getActiveCrystal();
	vec3 target = getAtom()->getAbsolutePosition();

	AtomGroupPtr chosen;
	chosen = crystal->getAtomsInBox(target, 5, 5, 5, true);
	_candidates.clear();

	for (int i = 0; i < chosen->atomCount(); i++)
	{
		if (chosen->atom(i)->getElementSymbol() == "H")
		{
			continue;
		}
		
		_candidates.push_back(chosen->atom(i));
	}
	
	if (_candidates.size() <= 3)
	{
		_disabled = true;
		return;
	}

	propagateChange();
}

AtomGroupPtr Sponge::correlGroup()
{
	AtomGroupPtr correl = AtomGroupPtr(new AtomGroup());
	correl->addAtom(getAtom());
	return correl;
}

void Sponge::addRestraintsToStrategy(RefinementStrategyPtr strategy)
{
	_anys.clear();
	for (Restraint::iterator it = _restraints.begin();
	      it != _restraints.end(); it++)
	{
		AnyPtr any = AnyPtr(new Any(&it->second));
		strategy->addParameter(&*any, Any::get, Any::set, 0.2, 0.01,
		shortDesc() + "_" + i_to_str(_anys.size()));
		_anys.push_back(any);
	}
}

std::string Sponge::shortDesc()
{
	return ("Sponge_" + getAtom()->shortDesc() + "_" +
	        i_to_str(_restraints.size()) + "r");
}

void Sponge::initialConnections()
{
	double bestScore = 0;
	AtomGroupPtr best = _close;
	
	for (int i = 0; i < 10; i++)
	{
		randomConnections();
		double score = AtomGroup::scoreWithMapGeneral(&_workspace, false);
		
		if (score < bestScore)
		{
			best = _close;
			bestScore = score;
		}
	}

	_close = best;
	std::cout << "Score now " << bestScore << ", score was "
	<< _preScore << std::endl;
	setupRestraints();
	singleRefine();
}

void Sponge::setupRestraints()
{
	for (int i = 0; i < _close->atomCount(); i++)
	{
		vec3 init = _close->atom(i)->getInitialPosition();
		vec3_subtract_from_vec3(&init, _abs);
		
		double length = vec3_length(init);
		_restraints[_close->atom(i)] = length;
	}
}

void Sponge::randomConnections()
{
	_restraints.clear();
	_close = AtomGroupPtr(new AtomGroup());
	std::random_shuffle(_candidates.begin(), _candidates.end());

	for (int i = 0; i < _candidates.size(); i++)
	{
		if (!_candidates[i]->getModel()->hasExplicitPositions())
		{
			continue;
		}
		
		if (_candidates[i] == getAtom())
		{
			continue;
		}
		
		if (_candidates[i]->getModel()->isSponge())
		{
			SpongePtr sp = ToSpongePtr(_candidates[i]->getModel());
			if (this != &*sp && sp->hasConnectedAtom(getAtom()))
			{
				continue;
			}
		}
		
		_close->addAtom(_candidates[i]);
		
		if (_close->atomCount() >= 4)
		{
			break;
		}
	}
	
	setupRestraints();
	singleRefine(true);
}

Sponge::Sponge(AtomPtr water) : Novalent(water)
{
	_disabled = false;
	_close = AtomGroupPtr(new AtomGroup());
	
	findConnections();
	setup_space(&_workspace);
	_workspace.crystal = Options::getActiveCrystal();
	_workspace.selectAtoms = correlGroup();
	
	/* old model */
	_preScore = AtomGroup::scoreWithMapGeneral(&_workspace, false);
	
	Novalent::getManyPositions();
	getFinalPositions();
	propagateChange();
}

void Sponge::singleRefine(bool others)
{
	propagateChange();
	Novalent::getManyPositions();
	propagateChange();

	getNetwork()->calculateSingle(shared_from_this(), false);
	
	if (others)
	{
		getNetwork()->calculateSingle(shared_from_this(), true);
	}
}

size_t Sponge::connectedCount()
{
	return _close->atomCount();
}

AtomPtr Sponge::getConnectedAtom(int i)
{
	return _close->atom(i);
}

Sponge::Sponge() : Novalent()
{
	_disabled = false;

}

Sponge::~Sponge()
{

}

void Sponge::copyActiveToFinalPos()
{
	std::lock_guard<std::mutex> lock(guiLock);

	_finalPositions[_n] = _storedSamples[_n].start;
}

std::vector<BondSample> *Sponge::getManyPositions(void *object)
{
	if (!_changedSamples)
	{
		return &_storedSamples;
	}

	/* should switch _changedSamples to false */
	getNetwork()->calculateSingle(shared_from_this());
	
	return &_storedSamples;
}

WaterNetworkPtr Sponge::getNetwork()
{
	AtomPtr a = getAtom();
	if (a->getMolecule() && a->getMolecule()->isWaterNetwork())
	{
		return ToWaterNetworkPtr(a->getMolecule());
	}

	return WaterNetworkPtr();
}

bool Sponge::hasConnectedAtom(AtomPtr atom)
{
	return _close->hasAtom(atom);
}

void Sponge::addAnglesToMap(std::map<AtomPtr,Restraint> &orig, 
                            std::map<AtomPtr,Restraint> &r)
{
	AtomPtr me = getAtom();
	double mean_angle = deg2rad(109.5);
	double stdev_angle = deg2rad(60);

	for (size_t i = 0; i < _candidates.size(); i++)
	{
		AtomPtr a = _candidates[i];
		double agree = orig[me][a];
		std::string symbol = a->getElementSymbol();
		
		if (symbol != "O" && symbol != "N")
		{
			continue;
		}


		vec3 apos = a->getInitialPosition();
		vec3_subtract_from_vec3(&apos, _abs);

		for (size_t j = 0; j < i; j++)
		{
			AtomPtr b = _candidates[j];
			symbol = b->getElementSymbol();

			if (symbol != "O" && symbol != "N")
			{
				continue;
			}

			double bgree = orig[me][b];

			vec3 bpos = b->getInitialPosition();
			vec3_subtract_from_vec3(&bpos, _abs);

			double angle = vec3_angle_with_vec3(apos, bpos);
			angle -= mean_angle;
			angle *= stdev_angle;

			double angree = exp(-(angle * angle));
			angree *= agree * bgree;

			angree = pow(angree, 1 / 3.);
			
//			std::cout << agree << " " << bgree << " " << angree << std::endl;
			
			if ((bgree < 0 && agree < 0) || angree != angree)
			{
				continue;
			}
			
			r[me][a] += angree;
			r[me][b] += angree;
			r[a][me] += angree;
			r[b][me] += angree;
			
			if (a->isHeteroAtom() || b->isHeteroAtom())
			{
				r[a][b] += angree;
				r[b][a] += angree;
			}
		}
	}
}

void Sponge::addLengthsToMap(std::map<AtomPtr,Restraint> &r)
{
	AtomPtr me = getAtom();

	double mean_length = 3.0;
	double stdev_length = 0.4;
	double penalty_length = 2.6;
	
	for (size_t i = 0; i < _candidates.size(); i++)
	{
		AtomPtr a = _candidates[i];
		
		if (me == a)
		{
			continue;
		}
		
		if (r.count(me) && r[me].count(a))
		{
			continue;
		}
		
		std::string symbol = a->getElementSymbol();
		
		if (symbol != "O" && symbol != "N")
		{
			continue;
		}

		vec3 apos = a->getInitialPosition();
		vec3_subtract_from_vec3(&apos, _abs);

		double length = vec3_length(apos);
		double orig = length;
		length -= mean_length;
		length *= stdev_length;
		
		double agree = exp(-length * length);
		
		if (orig < penalty_length)
		{
			double squash = orig * stdev_length;
			agree = exp(-squash * squash);
		}
		
		if (agree != agree)
		{
			continue;
		}
		
		r[me][a] = agree;
		r[a][me] = agree;
	}
}
