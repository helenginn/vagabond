// vagabond
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

#include <climits>
#include "ConfAxis.h"
#include "Bond.h"
#include "Atom.h"
#include "SVDBond.h"

ConfAxis::ConfAxis()
{
	_min = INT_MAX;
	_max = -INT_MAX;

}

void ConfAxis::getAxis(SVDBond *svd, int axis)
{
	if (axis >= svd->bondCount())
	{
		return;
	}

	for (size_t i = 0; i < svd->bondCount(); i++)
	{
		BondPtr b = svd->bond(i);
		AtomPtr a = b->getAtom();
		int resi = a->getResidueNum();
		
		double val = svd->svdValue(axis, i);
		val *= 2;

		setPhiDeviation(resi, val);
		setPsiDeviation(resi, -val);
	}

}
