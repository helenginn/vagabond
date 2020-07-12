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

#include "SymMate.h"
#include "Crystal.h"
#include "Anchor.h"
#include "Absolute.h"

SymMate::SymMate(CrystalPtr cryst)
{
	_cryst = cryst;
	_centre = cryst->centroid();
	_rot = make_mat3x3();
	_offset1 = empty_vec3();
	_offset2 = empty_vec3();

	mat3x3 recip = _cryst->getReal2Frac();
	std::cout << "Centroid: " << vec3_desc(_centre) <<  std::endl;
	mat3x3_mult_vec(recip, &_centre);
}

void closestToCentre(vec3 *v)
{
	while (v->x < -0.5) v->x += 1;
	while (v->x >= 0.5) v->x -= 1;

	while (v->y < -0.5) v->y += 1;
	while (v->y >= 0.5) v->y -= 1;

	while (v->z < -0.5) v->z += 1;
	while (v->z >= 0.5) v->z -= 1;
}

void rotateVec3(vec3 *v, mat3x3 rot, vec3 offset1, vec3 offset2)
{
	vec3 trial = *v;
	vec3_subtract_from_vec3(&trial, offset1);
	mat3x3_mult_vec(rot, &trial);
	vec3_add_to_vec3(&trial, offset1);
	vec3_subtract_from_vec3(&trial, offset2);
	*v = trial;
}

void SymMate::findSymop(vec3 target)
{
	std::cout << "Target: " << vec3_desc(target) << std::endl;
	mat3x3 recip = _cryst->getReal2Frac();
	mat3x3_mult_vec(recip, &target);

	CSym::CCP4SPG *spg = _cryst->getSpaceGroup();
	double min = FLT_MAX;
	int chosen = -1;
	vec3 first_offset = empty_vec3();
	vec3 second_offset = empty_vec3();
	vec3 empty = empty_vec3();
	
	for (int i = 0; i < spg->nsymop; i++)
	{
		mat3x3 rot = mat3x3_from_ccp4(spg->symop[i]);

		vec3 trial = _centre;
		vec3_subtract_from_vec3(&trial, target);

		vec3 original = trial;
		closestToCentre(&trial);
		vec3 running_offset = vec3_subtract_vec3(original, trial);

		mat3x3_mult_vec(rot, &trial);
		vec3_add_to_vec3(&trial, running_offset);

		original = trial;
		closestToCentre(&trial);

		double length = vec3_length(trial);
		
		if (length < min)
		{
			min = length;
			chosen = i;
			first_offset = running_offset;
			second_offset = vec3_subtract_vec3(original, trial);
			std::cout << "New length: " << min << std::endl;
		}
	}

	if (chosen < 0)
	{
		std::cout << "Finding closest symop failed." << std::endl;
		return;
	}

	_rot = mat3x3_from_ccp4(spg->symop[chosen]);
	_offset1 = first_offset;
	_offset2 = second_offset;
	std::cout << "First Offset: " << vec3_desc(_offset1) <<  std::endl;
	std::cout << "Rotation: " << std::endl;
	std::cout << mat3x3_desc(_rot) << std::endl;
	std::cout << "Second Offset: " << vec3_desc(_offset2) <<  std::endl;
}

void SymMate::applySymops(AtomGroupPtr group)
{
	vec3 oldcentre = group->centroid();
	mat3x3 recip = _cryst->getReal2Frac();
	mat3x3 real = _cryst->getHKL2Frac();

	vec3 copy = oldcentre;
	mat3x3_mult_vec(recip, &copy);
	std::cout << "Copy: " << vec3_desc(copy) << std::endl;
	rotateVec3(&copy, _rot, _offset1, _offset2);
	std::cout << "After: " << vec3_desc(copy) << std::endl;
	mat3x3_mult_vec(real, &copy);

	for (int i = 0; i < group->atomCount(); i++)
	{
		ModelPtr model = group->atom(i)->getModel();
		
		if (model->isBond())
		{
			continue;
		}
		
		if (model->isAbsolute())
		{
			AbsolutePtr abs = ToAbsolutePtr(model);
			vec3 pos = abs->getAbsolutePosition();
			vec3 init = group->atom(i)->getPDBPosition();
			
			mat3x3_mult_vec(recip, &pos);
			rotateVec3(&pos, _rot, _offset1, _offset2);
			mat3x3_mult_vec(real, &pos);

			mat3x3_mult_vec(recip, &init);
			rotateVec3(&init, _rot, _offset1, _offset2);
			mat3x3_mult_vec(real, &init);
			
			mat3x3 tensor = abs->getRealSpaceTensor();
			tensor = mat3x3_mult_mat3x3(_rot, tensor);
			abs->setTensor(tensor);

			abs->setPosition(pos);
			group->atom(i)->setPDBPosition(init);
		}

		/* not implementing anchor */
	}

	vec3 newcentre = group->centroid();
	std::cout << "Centre from " << vec3_desc(oldcentre) << " to " 
	<< vec3_desc(newcentre) << std::endl;
	std::cout << std::endl;
}
