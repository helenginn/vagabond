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

#include "SymAtom.h"
#include "Options.h"

SymAtom::SymAtom(Atom &parent) : Atom(parent)
{
	_parent = parent.shared_from_this();
	_symop = 0;
	_crystal = Options::getActiveCrystal();
	_ucShift = make_vec3(1, 1, 1);
	mat3x3 hkl2real = _crystal->getFrac2Real();
	mat3x3_mult_vec(hkl2real, &_ucShift);
}

SymAtom::~SymAtom()
{

}

vec3 SymAtom::transformVec(vec3 pos)
{
	vec3 v = _parent->getSymRelatedPosition(_symop, pos);
	vec3 ucDims = _ucShift;
	ucDims.x *= _shift.x;
	ucDims.y *= _shift.y;
	ucDims.z *= _shift.z;
	vec3_add_to_vec3(&v, ucDims);
	return v;
}

vec3 SymAtom::getAbsolutePosition()
{
	vec3 pos = _parent->getAbsolutePosition();
	return transformVec(pos);
}

vec3 SymAtom::getInitialPosition()
{
	vec3 pos = _parent->getInitialPosition();
	return transformVec(pos);
}
