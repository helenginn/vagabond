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

#include "KeyPoints.h"
#include "Param.h"
#include "Timer.h"
#include "Atom.h"
#include "RefinementNelderMead.h"
#include "Options.h"
#include "FlexGlobal.h"
#include "Anchor.h"
#include "Bond.h"
#include "Polymer.h"
#include <iostream>

KeyPoints::KeyPoints()
{

}

void KeyPoints::setPolymer(PolymerPtr polymer)
{
	_polymer = polymer;

	for (int i = 0; i < polymer->atomCount(); i++)
	{
		AtomPtr a = polymer->atom(i);
		
		if (a->getAtomName() != "CA")
		{
			continue;
		}
		
		ModelPtr mFirst = a->getModel();
		
		if (!mFirst->isBond())
		{
			continue;
		}
		
		BondPtr first = ToBondPtr(mFirst);
		first->setKeyPoints(shared_from_this());

		BondPtr second = first->downstreamBond(0, 0);
		
		if (!second)
		{
			continue;
		}

		second->setKeyPoints(shared_from_this());
	}
	
	int begin = _polymer->monomerBegin();
	int end = _polymer->monomerEnd();
	int step = (end - begin) / 10 + 1;
	
	for (int i = begin - step; i <= end + step; i += step)
	{
		WayPoint way;
		Param::setValue(&way.res, i);
		Param::setValue(&way.phi, 0);
		Param::setValue(&way.psi, 0);
		Param::setValue(&way.kick, 0);
		Param::setValue(&way.whack, 0);
		_points.push_back(way);
	}
}

double KeyPoints::getPhiContribution(BondPtr bond)
{
	return getContribution(bond, WayPointPhi);
}

double KeyPoints::getPsiContribution(BondPtr bond)
{
	return getContribution(bond, WayPointPsi);
}

double KeyPoints::getKickContribution(BondPtr bond)
{
	return getContribution(bond, WayPointKick);
}

double KeyPoints::getWhackContribution(BondPtr bond)
{
	return getContribution(bond, WayPointWhack);

}

double KeyPoints::getContribution(BondPtr bond, WayPointType type)
{
	int res = bond->getMinor()->getResidueNum();
	int wp = _points.size() - 2;
	
	for (int i = 0; i < _points.size() - 1; i++)
	{
		double before = Param::getValue(&_points[i].res);
		double after = Param::getValue(&_points[i + 1].res);
		
		if (after > res && before <= res)
		{
			wp = i;
			break;
		}
	}
	
	if (wp < 0)
	{
		wp = 0;
	}
	if (wp >= _points.size() - 1)
	{
		wp = _points.size() - 2;
	}
	
	double before = Param::getValue(&_points[wp].res);
	double after = Param::getValue(&_points[wp + 1].res);

	double bValue = 0; double aValue = 0;

	if (type == WayPointPhi)
	{
		bValue = _points[wp].phi.value();
		aValue = _points[wp + 1].phi.value();
	}
	else if (type == WayPointPsi)
	{
		bValue = _points[wp].psi.value();
		aValue = _points[wp + 1].psi.value();
	}
	else if (type == WayPointKick)
	{
		bValue = _points[wp].kick.value();
		aValue = _points[wp + 1].kick.value();
	}
	else if (type == WayPointWhack)
	{
		bValue = _points[wp].whack.value();
		aValue = _points[wp + 1].whack.value();
	}
	
	double prop = (after - res) / (after - before);
	double base = bValue;
	base += prop * (aValue - bValue);
	if (base != base)
	{
		base = 0;
	}

	return base;
}

bool KeyPoints::refineKeyPoints()
{
	CrystalPtr crystal = Options::getActiveCrystal();
	AtomGroupPtr backbone = _polymer->getAllBackbone();
	_global = FlexGlobal();
	_global.setAtomGroup(backbone);
	_global.setCrystal(crystal);
	
	double step = deg2rad(45);
	double tol = deg2rad(1);
	bool changed = false;

	std::cout << "Refining flexibility keypoints." << std::endl;

	{
		NelderMeadPtr nelder = NelderMeadPtr(new RefinementNelderMead());
		nelder->setCycles(30);
		nelder->setVerbose(true);	
		nelder->setSilent(true);

		Timer timer;

		for (int i = 0; i < _points.size(); i++)
		{
			nelder->addParameter(&_points[i].phi, Param::getValue, 
			                     Param::setValue, step, tol);
			nelder->addParameter(&_points[i].psi, Param::getValue, 
			                     Param::setValue, step, tol);
		}

		nelder->setEvaluationFunction(score, this);
		nelder->refine();
		nelder->reportResult();
		timer.quickReport();
		changed |= nelder->didChange();
	}

	std::cout << std::endl;
	return changed;

	{
		NelderMeadPtr nelder = NelderMeadPtr(new RefinementNelderMead());
		nelder->setCycles(60);
		nelder->setVerbose(true);	
		nelder->setSilent(true);

		Timer timer;
		double fac = (double)_polymer->monomerCount() / 50.;

		for (int i = 0; i < _points.size(); i++)
		{
			nelder->addParameter(&_points[i].kick, Param::getValue, 
			                     Param::setValue, 0.005 / fac, 0.0001 / fac);
		}

		nelder->setEvaluationFunction(score, this);
		nelder->refine();
		nelder->reportResult();
		timer.quickReport();
		changed |= nelder->didChange();
	}
	
	std::cout << std::endl;
	return changed;
}

double KeyPoints::score(void *object)
{
	KeyPoints *me = static_cast<KeyPoints *>(object);
	PolymerPtr pol = me->_polymer;
	pol->getAnchorModel()->forceRefresh();

	return me->_global.score(&me->_global);
}

void KeyPoints::addProperties()
{
	_tmpWays.clear();
	_tmpKicks.clear();

	for (int i = 0; i < _points.size(); i++)
	{
		vec3 pt = empty_vec3();
		pt.x = _points[i].res.value();
		pt.y = _points[i].phi.value();
		pt.z = _points[i].psi.value();
		_tmpWays.push_back(pt);
		_tmpKicks.push_back(make_vec3(_points[i].kick.value(),
		                              _points[i].whack.value(), 0));
	}

	addVec3ArrayProperty("waypoints", &_tmpWays);
	addVec3ArrayProperty("waykicks", &_tmpKicks);
	addReference("polymer", _polymer);

}

void KeyPoints::addObject(ParserPtr object, std::string category)
{

}

void KeyPoints::linkReference(BaseParserPtr object, std::string category)
{
	if (category == "polymer")
	{
		PolymerPtr polymer = ToPolymerPtr(object);
		_polymer = polymer;
	}
}

void KeyPoints::postParseTidy()
{
	_points.clear();
	
	for (int i = 0; i < _tmpWays.size(); i++)
	{
		WayPoint pt;
		pt.res.set_value(_tmpWays[i].x);
		pt.phi.set_value(_tmpWays[i].y);
		pt.psi.set_value(_tmpWays[i].z);
		pt.kick.set_value(0);
		
		if (_tmpKicks.size() > i)
		{
			pt.kick.set_value(_tmpKicks[i].x);
			pt.whack.set_value(_tmpKicks[i].y);
		}

		_points.push_back(pt);
	}
}
