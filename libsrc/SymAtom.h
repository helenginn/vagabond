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

#ifndef __vagabond__symatom__
#define __vagabond__symatom__

#include "shared_ptrs.h"
#include "Atom.h"

class SymAtom : public Atom
{
public:
	SymAtomPtr shared_from_this()
	{
		return ToSymAtomPtr(Parser::shared_from_this());
	}
	
	SymAtom(Atom &parent);
	
	virtual std::string shortDesc()
	{
		return Atom::shortDesc() + "sym";
	}
	
	void setSymop(int symop)
	{
		_symop = symop;
	}
	
	void setUnitCellShift(vec3 shift)
	{
		_shift = shift;
	}

	virtual vec3 getAbsolutePosition();
	virtual vec3 getInitialPosition();

	virtual ~SymAtom();
private:
	vec3 transformVec(vec3 pos);
	AtomPtr _parent;
	int _symop;
	vec3 _shift;
	vec3 _ucShift;
	CrystalPtr _crystal;
};

#endif
