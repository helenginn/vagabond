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

#include "Chelate.h"
#include "Atom.h"
#include "Element.h"
#include "Shouter.h"
#include "Absolute.h"

Chelate::Chelate()
{
	_realSpaceTensor = make_mat3x3();
	mat3x3_mult_scalar(&_realSpaceTensor, -1);
}

mat3x3 Chelate::getRealSpaceTensor()
{
	if (mat3x3_trace(_realSpaceTensor) < 0)
	{
		return _chelateAbs->getRealSpaceTensor();
	}

	return _realSpaceTensor;
}

void Chelate::setChelatedAtom(AtomPtr atom)
{
	if (!atom->getModel()->isAbsolute())
	{
		shout_at_helen("Chelated Atom is not Absolute!");
	}

	_atom = atom;
	_chelateAbs = ToAbsolutePtr(_atom->getModel());
	_atom->setModel(shared_from_this());
	refreshPositions();
}

void Chelate::addChelatingAtom(AtomPtr atom)
{
	if (!atom || atom == _atom)
	{
		return;
	}

	_chelating.push_back(atom);
	refreshPositions();
}

void Chelate::refreshPositions()
{
	if (!_atom || _chelating.size() == 0)
	{
		return;
	}
	
	mat3x3 empty = make_mat3x3();
	mat3x3_mult_scalar(&empty, 0);
	
	for (int i = 0; i < _chelating.size(); i++)
	{
		mat3x3 mat = _chelating[i]->getModel()->getRealSpaceTensor();
		mat3x3_add_mat3x3(&empty, mat);
	}
	
	double divide = _chelating.size();
	mat3x3_mult_scalar(&empty, 1 / divide);	
	
	empty.vals[0] += b2var(_bufferB);
	empty.vals[4] += b2var(_bufferB);
	empty.vals[8] += b2var(_bufferB);
	
	_realSpaceTensor = empty;
	_chelateAbs->setTensor(_realSpaceTensor);
}

double Chelate::getMeanSquareDeviation()
{
	if (_chelating.size())
	{
		double trace = mat3x3_trace(_realSpaceTensor);
		trace = var2b(trace);
		
		if (trace > 0)
		{
			return trace;
		}
	}

	return _chelateAbs->getMeanSquareDeviation();
}

FFTPtr Chelate::makeDistribution()
{
	return _chelateAbs->makeDistribution();
}

vec3 Chelate::getAbsolutePosition()
{
	return _chelateAbs->getAbsolutePosition();
}

std::string Chelate::shortDesc()
{
	std::string str = "Ch" + _atom->getElement()->getSymbol();
	str += i_to_str(_atom->getAtomNum());
	str += "_" + i_to_str(_chelating.size());
	
	return str;
}
