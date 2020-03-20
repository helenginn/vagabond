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
#include "Any.h"
#include "Options.h"
#include "ExplicitModel.h"
#include "WaterNetwork.h"

void Sponge::findConnections()
{
	CrystalPtr crystal = Options::getActiveCrystal();
	std::vector<AtomPtr> chosen = crystal->getCloseAtoms(getAtom(), 5, false);
	_candidates.clear();

	for (int i = 0; i < chosen.size(); i++)
	{
		if (chosen[i]->getElementSymbol() == "H")
		{
			continue;
		}
		
		_candidates.push_back(chosen[i]);
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
	for (Restraint::iterator it = _restraints.begin();
	      it != _restraints.end(); it++)
	{
		AnyPtr any = AnyPtr(new Any(&it->second));
		strategy->addParameter(&*any, Any::get, Any::set, 0.2, 0.01);
	}
}

std::string Sponge::shortDesc()
{
	return ("Sponge_" + getAtom()->shortDesc() + 
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
	getNetwork()->recalculate();
	
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

