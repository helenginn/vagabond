// Vagabond
// Copyright (C) 2017-2018 Helen Ginn
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

#ifndef __Vagabond__Balance__h
#define __Vagabond__Balance__h

#include "shared_ptrs.h"
#include "RefinementStrategy.h"

/**
 * \class Balance
 * \brief Looks after a balancing of occupancy parameters which all
 * have to add up to 1. */

class Balance
{
public:
	Balance(BondPtr bond);

	double addParamsToStrategy(RefinementStrategyPtr strategy);

	bool isFromBond(BondPtr bond)
	{
		return (_bond == bond);
	}
	void adjustment();
private:
	std::vector<ParamBandPtr> _bands;
	
	void *_evalObj;
	Getter _evalFunc;
	BondPtr _bond;
};

#endif
