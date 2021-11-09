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

#ifndef __vagabond__ConfAxis__
#define __vagabond__ConfAxis__

#include "shared_ptrs.h"
#include <map>
#include <iostream>

class ConfAxis
{
public:
	ConfAxis();
	
	void setWhackDeviation(int resi, double val)
	{
		_residueWhacks[resi] = val;
	}
	
	void setTorsionDeviation(int resi, double val)
	{
		_residueTorsions[resi] = val;
	}

	double getTorsionDeviationForResidue(int i)
	{
		if (_residueTorsions.count(i))
		{
			return _residueTorsions[i];
		}
		
		return 0.0;
	}

	double getWhackDeviationForResidue(int i)
	{
		if (_residueWhacks.count(i))
		{
			return _residueWhacks[i];
		}
		
		return 0.0;
	}
	
private:
	std::map<int, double> _residueTorsions;
	std::map<int, double> _residueWhacks;

};

#endif
