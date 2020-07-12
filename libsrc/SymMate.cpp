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
#include "FFT.h"
#include "Anchor.h"
#include "Absolute.h"

SymMate::SymMate(CrystalPtr cryst)
{
	_cryst = cryst;
	_centre = cryst->centroid();
	_rot = make_mat3x3();
	_offset = empty_vec3();

	mat3x3 recip = _cryst->getReal2Frac();
	mat3x3_mult_vec(recip, &_centre);
	std::cout << "Centroid: " << vec3_desc(_centre) <<  std::endl;
}

void SymMate::findSymop(vec3 target)
{
	mat3x3 recip = _cryst->getReal2Frac();
	mat3x3_mult_vec(recip, &target);

	CSym::CCP4SPG *spg = _cryst->getSpaceGroup();
	double min = FLT_MAX;
	int chosen = -1;
	vec3 offset = empty_vec3();
	
	for (int i = 0; i < spg->nsymop; i++)
	{
		vec3 trial = _centre;
		mat3x3 rot = mat3x3_from_ccp4(spg->symop[i]);
		mat3x3_mult_vec(rot, &trial);
		
		vec3_subtract_from_vec3(&trial, target);
		vec3 original = trial;
		VagFFT::collapseFrac(&trial.x, &trial.y, &trial.z);

		double length = vec3_length(trial);
		
		if (length < min)
		{
			min = length;
			chosen = i;
			offset = vec3_subtract_vec3(original, trial);
		}
	}

	if (chosen < 0)
	{
		std::cout << "Finding closest symop failed." << std::endl;
		return;
	}

	_rot = mat3x3_from_ccp4(spg->symop[chosen]);
	_offset = offset;
	std::cout << "Rotation: " << std::endl;
	std::cout << mat3x3_desc(_rot) << std::endl;
	std::cout << "Offset: " << vec3_desc(_offset) <<  std::endl;
}

void SymMate::applySymops(AtomGroupPtr group)
{
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
			
			mat3x3_mult_vec(_rot, &pos);
			vec3_subtract_from_vec3(&pos, _offset);
			
			mat3x3 tensor = abs->getRealSpaceTensor();
			tensor = mat3x3_mult_mat3x3(_rot, tensor);
			abs->setTensor(tensor);

			abs->setPosition(pos);
		}

		if (model->isAnchor())
		{
			AnchorPtr anch = ToAnchorPtr(model);
			anch->applyRotation(_rot);
			anch->applyOffset(_offset);
		}
	}

}
