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

#ifndef __vagabond__Superpose__
#define __vagabond__Superpose__

#include <hcsrc/mat3x3.h>
#include <hcsrc/vec3.h>
#include "shared_ptrs.h"
#include "ExplicitModel.h"
#include <map>
#include <vector>

class Superpose
{
public:
	Superpose();

	void savePositions();
	void calculateDeviations();
	
	void applyDeviations(std::vector<BondSample> &samples);
	
	void setAtoms(AtomList full)
	{
		_full = full;
	}
	
	size_t atomCount()
	{
		return _atoms.size();
	}
private:
	std::vector<AtomPtr> _atoms;
	std::vector<AtomPtr> _full;
	std::map<AtomPtr, std::vector<vec3> > _saved;
	
	std::vector<vec3> _remove;
	std::vector<vec3> _add;
	std::vector<mat3x3> _rotations;

	size_t _sampleSize;
};

#endif
