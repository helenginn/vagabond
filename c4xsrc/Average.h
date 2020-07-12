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

#ifndef __cluster4x__Average__
#define __cluster4x__Average__

#include "MtzFFTPtr.h"
#include <vector>

class Group;

class Average
{
public:
	Average(Group *group);

	virtual ~Average()
	{

	}

	virtual void calculate() = 0;
	virtual void findIntercorrelations(Group *other, double **svd);
protected:
	std::vector<MtzFFTPtr> _mtzs;
	bool _symmetric;
private:
	virtual double findCorrelation(MtzFFTPtr one, MtzFFTPtr two) = 0;
	Group *_group;
};

#endif
