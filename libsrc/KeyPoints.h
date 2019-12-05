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

#ifndef __Vagabond__KeyPoints__h
#define __Vagabond__KeyPoints__h

#include <vector>
#include "shared_ptrs.h"
#include "FlexGlobal.h"
#include "Param.h"

typedef struct
{
	Param res;
	Param phi;
	Param psi;
} WayPoint;

/** 
 * \class KeyPoints
 * \brief specify key points for a given monomer to derive phi/psi
 * modifications for a Bond */

class KeyPoints : public boost::enable_shared_from_this<KeyPoints>
{
public:
	KeyPoints();
	
	void setPolymer(PolymerPtr polymer);

	double getPhiContribution(BondPtr bond);
	double getPsiContribution(BondPtr bond);
	
	bool refineKeyPoints();
private:
	static double score(void *object);
	std::vector<WayPoint> _points;

	PolymerPtr _polymer;
	FlexGlobal _global;

};

#endif
