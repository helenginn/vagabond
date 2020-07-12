// cluster4x
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

#include "AveCAlpha.h"
#include "MtzFile.h"
#include "MtzFFT.h"
#include <libsrc/maths.h>
#include "QuickAtoms.h"
#include <iostream>

AveCAlpha::AveCAlpha(Group *group) : Average(group)
{

}

void AveCAlpha::calculate()
{
	_quick = new QuickAtoms(NULL);

	for (size_t i = 0; i < _mtzs.size(); i++)
	{
		if (_mtzs[i]->getMtzFile()->isDead())
		{
			continue;
		}

		QuickAtoms *quick = _mtzs[i]->getMtzFile()->getQuickAtoms();
		quick->fetchAtoms();
		_quick->addAtomsFrom(quick);
	}

	_quick->divideThrough();
}

double AveCAlpha::findCorrelation(MtzFFTPtr one, MtzFFTPtr two)
{
	QuickAtoms *qOne = one->getMtzFile()->getQuickAtoms();
	QuickAtoms *qTwo = two->getMtzFile()->getQuickAtoms();
	
	double cc = 0;
	cc = QuickAtoms::compare(qOne, qTwo, _quick);
	return cc;
}

vec3 AveCAlpha::getCentre()
{
	if (!_quick)
	{
		return empty_vec3();
	}

	return _quick->getCentre();
}
