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

#include "Motion.h"
#include "Anisotropicator.h"

Motion::Motion()
{
	_trans = RefineMat3x3Ptr(new RefineMat3x3(this, NULL));

}

void Motion::translateStartPositions(std::vector<BondSample> &stored)
{
	mat3x3 translation = _trans->getMat3x3();
	Anisotropicator tropicator;
	tropicator.setTensor(translation);
	mat3x3 trans = tropicator.basis();

	vec3 sum_start = empty_vec3();
	vec3 sum_old = empty_vec3();

	for (int i = 0; i < stored.size(); i++)
	{
		vec3_add_to_vec3(&sum_start, stored[i].start);
		vec3_add_to_vec3(&sum_old, stored[i].old_start);
	}
	
	vec3_mult(&sum_start, 1 / (double)stored.size());
	vec3_mult(&sum_old, 1 / (double)stored.size());

	for (int i = 0; i < stored.size(); i++)
	{
		vec3 start = stored[i].start;
		vec3 diff = vec3_subtract_vec3(start, sum_start);
		mat3x3_mult_vec(trans, &diff);
		vec3 new_start = vec3_add_vec3(sum_start, diff);
		vec3 old_start = vec3_add_vec3(sum_old, diff);

		stored[i].start = new_start;
		stored[i].old_start = old_start;
	}
}


