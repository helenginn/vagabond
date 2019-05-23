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
#include <iomanip>
#include "Bond.h"
#include "Anchor.h"
#include "Shouter.h"
#include "mat3x3.h"
#include "vec3.h"

Whack::Whack()
{
	_kick = 0;
	_whack = 0;
	_valid = true;
	_enabled = true;
}

void Whack::setBond(BondPtr bond)
{
	_bond = bond;
	
	/* Check that the following bond is allowed to refine flexibility */
	if (_bond->downstreamBondGroupCount() && _bond->downstreamBondCount(0))
	{
		BondPtr child = _bond->downstreamBond(0, 0);
		
		if (!child->getRefineFlexibility())
		{
			_valid = false;
			return;
		}
	}
	else
	{
		_valid = false;
		return;
	}
	
	_bond->setWhack(shared_from_this());
	saveSamples();
	applyKick();
}

bool Whack::needsRefresh(std::vector<BondSample> &anchSamp)
{
	return (anchSamp.size() != _samples.size());
}

void Whack::saveSamples()
{
	_samples = *_bond->getManyPositions();
	BondPtr child = _bond->downstreamBond(0, 0);
	child->calculateMagicMat();
	child->correctTorsionAngles();
	std::vector<BondSample> chSamples = *child->getManyPositions();
	
	for (size_t i = 0; i < chSamples.size(); i++)
	{
		_samples[i].kickValue = chSamples[i].kickValue;
	}

}

void Whack::applyKick()
{
	if (!_bond)
	{
		shout_at_helen("Whack doesn't have bond set.");
	}
	
	if (!_valid)
	{
		return;
	}
	
	if (_bond->downstreamBondGroupCount() && _bond->downstreamBondCount(0))
	{
		double value = 1 * _whack + _kick;
		
		if (!_enabled)
		{
			value = 0;
		}

		BondPtr child = _bond->downstreamBond(0, 0);
		Bond::setKick(&*child, value);
//		child->propagateChange(-1);
	}
}

void Whack::setWhack(void *object, double whack)
{
	Whack *me = static_cast<Whack *>(object);
	me->_whack = whack;
	me->applyKick();
}

void Whack::addToAnchor(AnchorPtr anchor)
{
	if (_valid)
	{
		_anchor = anchor;
		anchor->addWhack(shared_from_this());
	}
}

void Whack::applyToAnchorSamples(std::vector<BondSample> &anchSamp)
{
	if (fabs(_whack) < 1e-6 || !_enabled)
	{
		return;
	}
	
	AnchorPtr anchor = _anchor.lock();
	
	AtomPtr anchAtom = anchor->getAtom();
	AtomPtr bondAtom = _bond->getAtom();
	
	double check = 0;

	for (int i = 0; i < anchSamp.size(); i++)
	{
		double kickvalue = _samples[i].kickValue;
		double mag = kickvalue * _whack;
		anchSamp[i].kickValue = mag;
		check += mag;
	}

	check /= (double)anchSamp.size();
	
	for (int i = 0; i < anchSamp.size(); i++)
	{
		mat3x3 bond_basis = _samples[i].basis;
		vec3 bond_axis = mat3x3_axis(bond_basis, 2);

		double mag = anchSamp[i].kickValue;
		mat3x3 rot = mat3x3_unit_vec_rotation(bond_axis, -mag + check);

		vec3 start = anchSamp[i].start;
		vec3 old_start = anchSamp[i].old_start;
		vec3 centre = _samples[i].start;
		
		anchSamp[i].start = rotate_round_bond(start, centre, rot);
		anchSamp[i].old_start = rotate_round_bond(old_start, centre, rot);
		
		mat3x3 basis = anchSamp[i].basis;
		anchSamp[i].basis = mat3x3_mult_mat3x3(rot, basis);
	}
}

void Whack::setKick(void *object, double kick)
{
	Whack *me = static_cast<Whack *>(object);
	me->_kick = kick;
	me->applyKick();
}

std::string Whack::getParserIdentifier()
{
	return "whack_" + _bond->shortDesc();
}

void Whack::addProperties()
{
	addDoubleProperty("whack", &_whack);
	addDoubleProperty("kick", &_kick);
	addReference("anchor", _anchor.lock());
	addReference("bond", _bond);
}

void Whack::linkReference(BaseParserPtr object, std::string category)
{
	if (category == "anchor")
	{
		_anchor = ToAnchorPtr(object);
	}
	else if (category == "bond")
	{
		_bond = ToBondPtr(object);
	}
}

void Whack::postParseTidy()
{

}
