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

#ifndef __vagabond__SymMate__
#define __vagabond__SymMate__

#include "shared_ptrs.h"
#include "vec3.h"
#include "mat3x3.h"

class SymMate
{
public:
	SymMate(CrystalPtr cryst);

	void findSymop(vec3 target);
	void applySymops(AtomGroupPtr group);
private:
	CrystalPtr _cryst;
	vec3 _realcentre;
	vec3 _centre;
	vec3 _target;

	mat3x3 _rot;
	vec3 _screw;
	vec3 _offset1;
	vec3 _offset2;
};


#endif
