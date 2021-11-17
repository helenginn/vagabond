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

class SVDBond;

class ConfAxis
{
public:
	ConfAxis();
	
	void getAxis(SVDBond *svd, int i);
	
	void setPhiDeviation(int resi, double val)
	{
		_residuePhis[resi] = val;
		if (_max < resi)
		{
			_max = resi;
		}
	}
	
	void setPsiDeviation(int resi, double val)
	{
		_residuePsis[resi] = val;
		if (_min > resi)
		{
			_min = resi;
		}
	}

	double getPsiDeviationForResidue(int i)
	{
		return _residuePsis[i];
	}

	double getPhiDeviationForResidue(int i)
	{
		return _residuePhis[i];
	}
	
	int residueBegin()
	{
		return _min;
	}

	int residueEnd()
	{
		return _max;
	}
	
	int count()
	{
		return _residuePsis.size();
	}
private:
	std::map<int, double> _residuePsis;
	std::map<int, double> _residuePhis;
	int _min;
	int _max;

};

#endif
