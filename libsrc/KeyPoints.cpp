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
	const int offset = 3;
	const int step = 6;
	
	for (int i = begin + 3; i < end; i += step)
	{
		WayPoint way;
		Param::setValue(&way.res, i);
		Param::setValue(&way.psi, 0);
		Param::setValue(&way.phi, 0);
		_points.push_back(way);
	}
}

double KeyPoints::getPhiContribution(BondPtr bond)
{
	int res = bond->getMinor()->getResidueNum();
	int wp = _points.size() - 1;
	
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
	
	double before = Param::getValue(&_points[wp].res);
	double bValue = Param::getValue(&_points[wp].phi);

	double after = Param::getValue(&_points[wp + 1].res);
	double aValue = Param::getValue(&_points[wp + 1].phi);
	
	double prop = (after - res) / (after - before);
	prop *= M_PI;
	double cosProp = cos(prop);
	
	cosProp += 1;
	double multiply = (bValue - aValue) / 2;
	cosProp *= multiply;
	cosProp += aValue;

	return cosProp;
}

double KeyPoints::getPsiContribution(BondPtr bond)
{
	return 0;
}

bool KeyPoints::refineKeyPoints()
{
	CrystalPtr crystal = Options::getActiveCrystal();
	AtomGroupPtr backbone = _polymer->getAllBackbone();
	_global = FlexGlobal();
	_global.setAtomGroup(backbone);
	_global.setCrystal(crystal);
	_global.matchElectronDensity();
	
	double step = deg2rad(45);
	double tol = deg2rad(1);

	NelderMeadPtr nelder = NelderMeadPtr(new RefinementNelderMead());
	nelder->setCycles(60);
	nelder->setVerbose(true);	
	nelder->setSilent(true);

	for (int i = 0; i < _points.size(); i++)
	{
		nelder->addParameter(&_points[i].phi, Param::getValue, 
		                     Param::setValue, step, tol);
	}

	nelder->setEvaluationFunction(score, this);
	nelder->refine();
	nelder->reportResult();
	
	return nelder->didChange();
}

double KeyPoints::score(void *object)
{
	KeyPoints *me = static_cast<KeyPoints *>(object);
	PolymerPtr pol = me->_polymer;
	pol->getAnchorModel()->forceRefresh();
	pol->propagateChange();

	return me->_global.score(&me->_global);
}
