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

#include "GhostBond.h"
#include "Atom.h"

GhostBond::GhostBond()
{

}

void GhostBond::setAtoms(AtomPtr major, AtomPtr minor)
{
	if (!major || !minor)
	{
		return;
	}

	_major = major;
	_minor = minor;
	
	minor->setGhostBond(shared_from_this());
}

std::string GhostBond::getParserIdentifier()
{
	return "ghost_" + i_to_str(getMinor()->getAtomNum()) +
	getMajor()->getAtomName() + "-" + getMinor()->getAtomName();
}

void GhostBond::addProperties()
{
	addReference("major", _major.lock());
	addReference("minor", _minor.lock());
}

void GhostBond::linkReference(ParserPtr object, std::string category)
{
	if (category == "major")
	{
		_major = ToAtomPtr(object);
	}
	else if (category == "minor")
	{
		_minor = ToAtomPtr(object);
	}
}

void GhostBond::postParseTidy()
{

}
