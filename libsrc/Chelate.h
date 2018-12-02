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

#ifndef __vagabond__Chelate__
#define __vagabond__Chelate__

#include "Model.h"
#include "mat3x3.h"

class Chelate : public Model
{
public:
	Chelate();

	virtual AtomPtr getAtom()
	{
		return _atom;
	}
	
	/* Get a new anisotropic tensor */
	virtual void refreshPositions();
	
	virtual bool hasExplicitPositions()
	{
		return false;
	}
	
	virtual std::string shortDesc();
	
	void setBufferB(double buf)
	{
		_bufferB = buf;
	}

	virtual mat3x3 getRealSpaceTensor();
	
	void setChelatedAtom(AtomPtr atom);
	void addChelatingAtom(AtomPtr atom);

	vec3 getAbsolutePosition();
	virtual double getMeanSquareDeviation();
	FFTPtr makeDistribution();
	
	virtual std::string getClassName()
	{
		return "Chelate";
	}
protected:

	virtual std::string getParserIdentifier()
	{
		return "chelate"; 
	}
private:
	double _bufferB;
	AtomPtr _atom;
	AbsolutePtr _chelateAbs;
	
	std::vector<AtomPtr> _chelating;

};

#endif
