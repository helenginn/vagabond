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

#ifndef __cluster4x__multistate__
#define __cluster4x__multistate__

#include "shared_ptrs.h"

class Multistate
{
public:
	Multistate(std::string filename);
	
	void ignoreAtomsExcept(std::string atom)
	{
		_atom = atom;
	}

	void process();
	
	size_t crystalCount()
	{
		return _crystals.size();
	}
	
	CrystalPtr crystal(int i)
	{
		return _crystals[i];
	}
private:
	std::string _filename;
	std::string _atom;
	std::string _cryst1;
	std::vector<CrystalPtr> _crystals;
	std::vector<std::string> _individuals;
};

#endif
