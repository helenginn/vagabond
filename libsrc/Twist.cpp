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

#include "Twist.h"
#include "Anchor.h"
#include "Bond.h"

Twist::Twist()
{
	_enabled = true;
	_valid = false;
	_twist = 0;
}

void Twist::setBond(BondPtr bond)
{
	if (bond->downstreamBondGroupCount() == 0)
	{
		return;
	}

	BondPtr child = bond->downstreamBond(0, 0);

	if (child->isFixed() || !child->isUsingTorsion())
	{
		return;
	}
	
	if (bond->hasTwist())
	{
		std::cout << "Already has twist" << std::endl;
		return;
	}
	
	_bond = bond;
	saveSamples();
	_valid = true;
	
	bond->setTwist(shared_from_this());
}

void Twist::addToAppliedModel(ExplicitModelPtr applied)
{
	if (_valid)
	{
		_applied = applied;
		applied->addTwist(shared_from_this());
	}
}

void Twist::saveSamples()
{
	_samples = *_bond->getManyPositions();
}


void Twist::applyToAnchorSamples(std::vector<BondSample> &anchSamp)
{
	if (fabs(_twist) < 1e-6 || !_enabled)
	{
		return;
	}

	for (int i = 0; i < anchSamp.size(); i++)
	{
		mat3x3 bond_basis = _samples[i].basis;
		vec3 bond_axis = mat3x3_axis(bond_basis, 2);

		mat3x3 rot = mat3x3_unit_vec_rotation(bond_axis, _twist);

		vec3 start = anchSamp[i].start;
		vec3 old_start = anchSamp[i].old_start;
		vec3 centre = _samples[i].start;
		
	
		anchSamp[i].start = rotate_round_bond(start, centre, rot);
		anchSamp[i].old_start = rotate_round_bond(old_start, centre, rot);

		mat3x3 basis = anchSamp[i].basis;
		anchSamp[i].basis = mat3x3_mult_mat3x3(rot, basis);
	}
}

void Twist::setTwist(void *object, double val)
{
	Twist *twist = static_cast<Twist *>(object);
	twist->_twist = val;
	twist->getAppliedModel()->propagateChange(20);
}

