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

#include "Superpose.h"
#include "Atom.h"
#include "Bond.h"
#include <hcsrc/Matrix.h>
#include <iomanip>

Superpose::Superpose()
{
	_sampleSize = 0;

}

void Superpose::savePositions()
{
	_saved.clear();
	_add.clear();

	for (size_t i = 0; i < _full.size(); i++)
	{
		if (!_full[i]->getModel()->hasExplicitPositions())
		{
			continue;
		}
		
		ExplicitModelPtr model = ToExplicitModelPtr(_full[i]->getModel());
		const std::vector<BondSample> *samples = model->getManyPositions();
		std::vector<vec3> positions;
		positions.reserve(samples->size());
		
		_sampleSize = samples->size();
		
		if (_sampleSize <= 5)
		{
			return;
		}
		
		for (size_t j = 0; j < samples->size(); j++)
		{
			vec3 p = samples->at(j).start;
			positions.push_back(p);
		}
		
		_saved[_full[i]] = positions;
		_atoms.push_back(_full[i]);
	}
	
	for (size_t i = 0; i < _sampleSize; i++)
	{
		vec3 avePos = empty_vec3();
		double count = 0;

		for (size_t j = 0; j < _atoms.size(); j++)
		{
			vec3 pos = _saved[_atoms[j]][i];
			avePos += pos;
			count++;
		}
		
		avePos /= count;
		_add.push_back(avePos);
	}
}

void Superpose::calculateDeviations()
{
	if (_atoms.size() == 0)
	{
		savePositions();
	}

	_remove.clear();
	_rotations.clear();

	for (size_t n = 0; n < _sampleSize; n++)
	{
		vec3 newAve = empty_vec3();
		double count = 0;
		bool stop = false;

		for (size_t i = 0; i < _atoms.size(); i++)
		{
			AtomPtr atom = _atoms[i];
			ExplicitModelPtr model = ToExplicitModelPtr(atom->getModel());
			std::vector<BondSample> *samples = model->getManyPositions();
			
			if (n >= samples->size())
			{
				stop = true;
				break;
			}

			vec3 update = samples->at(n).start;
			newAve += update;
			count++;
		}
		
		if (stop)
		{
			continue;
		}

		newAve /= count;
		_remove.push_back(newAve);

		HelenCore::SVD svd;
		setupSVD(&svd, 3, 3);
		for (size_t i = 0; i < _atoms.size(); i++)
		{
			AtomPtr atom = _atoms[i];
			ExplicitModelPtr model = ToExplicitModelPtr(atom->getModel());
			std::vector<BondSample> *samples = model->getManyPositions();

			vec3 update = samples->at(n).start;
			vec3 old = _saved[atom][n];
			update -= _add[n];
			old -= _remove[n];
			
			for (size_t k = 0; k < 3; k++)
			{
				for (size_t j = 0; j < 3; j++)
				{
					double add = *(&update.x + j) * *(&old.x + k);
					svd.u.ptrs[j][k] += add;
				}
			}
		}

		bool success = runSVD(&svd);
		
		mat3x3 rot = make_mat3x3();
		
		if (success)
		{
			for (size_t i = 0; i < 9; i++)
			{
				rot.vals[i] = 0;
			}

			for (size_t j = 0; j < 3; j++)
			{
				for (size_t i = 0; i < 3; i++)
				{
					for (size_t n = 0; n < 3; n++)
					{
						rot.vals[i+3*j] += svd.v.ptrs[i][n] * svd.u.ptrs[j][n];
					}
				}
			}
		}

		rot = mat3x3_transpose(rot);
		_rotations.push_back(rot);

		freeSVD(&svd);
	}
}

void Superpose::applyDeviations(std::vector<BondSample> &samples)
{
	if (samples.size() <= 5)
	{
		_atoms.clear();
		return;
	}

	for (size_t i = 0; i < samples.size() && i < _rotations.size(); i++)
	{
		samples[i].start -= _remove[i];
		mat3x3_mult_vec(_rotations[i], &samples[i].start);
		samples[i].start += _add[i];

		samples[i].old_start -= _remove[i];
		mat3x3_mult_vec(_rotations[i], &samples[i].old_start);
		samples[i].old_start += _add[i];
	}

}
