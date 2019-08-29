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

#ifndef __vagabond__Motion__
#define __vagabond__Motion__

#include "RefineMat3x3.h"
#include "ExplicitModel.h"

/**
 * \class Motion
 * \brief Whole-molecule motions which can be applied to multiple anchor
 * points if needed. 
 **/

class Motion
{
public:
	Motion();


private:
	void translateStartPositions(std::vector<BondSample> &stored);
	RefineMat3x3Ptr _trans;

};

#endif
