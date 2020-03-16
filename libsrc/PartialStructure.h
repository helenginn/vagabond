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

#ifndef __vagabond__partialstructure__
#define __vagabond__partialstructure__

#include "shared_ptrs.h"

/** \class PartialStructure
 *  \brief Will scale with absolute value and B factor applied to a set
 *   of partial structure factors which may either be input from the
 *   command line or created by Vagabond (e.g. solvent models).
 */

class PartialStructure
{
public:
	PartialStructure() {};
	
	void setStructure(VagFFTPtr refPart);
	void scalePartialStructure();

	static double getSolvScale(void *object)
	{
		return static_cast<PartialStructure *>(object)->_solvScale;
	}
	
	static void setSolvScale(void *object, double value)
	{
		static_cast<PartialStructure *>(object)->_solvScale = value;
	}

	static double getSolvBFac(void *object)
	{
		return static_cast<PartialStructure *>(object)->_solvBFac;
	}
	
	static void setSolvBFac(void *object, double value)
	{
		static_cast<PartialStructure *>(object)->_solvBFac = value;
	}

	void setData(DiffractionPtr data)
	{
		_data = data;
	}
	
	void setCrystal(CrystalPtr crystal)
	{
		_crystal = crystal;
	}
	
	CrystalPtr getCrystal()
	{
		return _crystal.lock();
	}

	virtual void reportScale();
protected:
	CrystalWkr _crystal;
	double _solvScale;
	double _solvBFac;

	void setPartialStructure(VagFFTPtr solvent)
	{
		_partial = solvent;
	}
private:
	DiffractionPtr _data;
	double scaleAndAddPartialScore();
	double scaleOrScore(bool score);
	static double scalePartialScore(void *object);
	VagFFTPtr _partial;
	
};

#endif
