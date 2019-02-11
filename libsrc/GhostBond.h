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

#ifndef __vagabond__ghostbond__
#define __vagabond__ghostbond__

#include "Parser.h"
#include "shared_ptrs.h"

class GhostBond : public Parser
{
public:
	GhostBond();

	GhostBondPtr shared_from_this()
	{
		return ToGhostBondPtr(Parser::shared_from_this());
	}

	void setAtoms(AtomPtr major, AtomPtr minor);
	
	AtomPtr getMajor()
	{
		return _major.lock();
	}
	
	AtomPtr getMinor()
	{
		return _minor.lock();
	}
	
	std::string getClassName()
	{
		return "Ghost";
	}
	
protected:
	virtual std::string getParserIdentifier();
	virtual void addProperties();
	virtual void linkReference(ParserPtr object, std::string category);
	virtual void postParseTidy();

private:
	AtomWkr _major;
	AtomWkr _minor;

};

#endif
