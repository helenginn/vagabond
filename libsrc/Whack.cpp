// Vagabond : bond-based macromolecular model refinement
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

#include "Whack.h"
#include "Bond.h"
#include "Anchor.h"
#include "Shouter.h"
#include "mat3x3.h"
#include "vec3.h"

Whack::Whack()
{
	_prop = 1;
	_whack = 0;
}

void Whack::setBond(BondPtr bond)
{
	_bond = bond;
	applyKick();
}

void Whack::applyKick()
{
	if (!_bond)
	{
		shout_at_helen("Whack doesn't have bond set.");
	}
	
	Bond::setKick(&*_bond, _whack * _prop);
	_samples = *_bond->getManyPositions();
}

void Whack::setWhack(void *object, double whack)
{
	Whack *me = static_cast<Whack *>(object);
	me->_whack = whack;
	me->applyKick();
}

void Whack::addToAnchor(AnchorPtr anchor)
{
	_anchor = anchor;
	anchor->addWhack(shared_from_this());
}

void Whack::applyToAnchorSamples(std::vector<BondSample> &anchSamp)
{
	for (int i = 0; i < anchSamp.size(); i++)
	{
		double kick = _samples[i].kickValue;
		double mag = kick * (1 - _prop) * _whack;
		mat3x3 bond_basis = _samples[i].basis;
		vec3 bond_axis = mat3x3_axis(bond_basis, 2);
		mat3x3 rot = mat3x3_unit_vec_rotation(bond_axis, -mag);

		vec3 start = anchSamp[i].start;
		anchSamp[i].start = mat3x3_mult_vec(rot, start);
		
		mat3x3 basis = anchSamp[i].basis;
		mat3x3 transpose = mat3x3_transpose(rot);

		anchSamp[i].basis = mat3x3_mult_mat3x3(transpose, basis);
	}
}

void Whack::setProportion(void *object, double prop)
{
	Whack *me = static_cast<Whack *>(object);
	me->_prop = prop;
	me->applyKick();
}
