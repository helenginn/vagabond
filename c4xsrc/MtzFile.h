// Clusterxxxx
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

#ifndef __clusterxxxx__MtzFile__
#define __clusterxxxx__MtzFile__

#include <libsrc/shared_ptrs.h>
#include <libsrc/vec3.h>
#include <string>
#include "GLObject.h"

class MtzFile
{
public:
	MtzFile(std::string filename);
	
	void setRefinementID(std::string id)
	{
		_refinementID = id;
	}
	
	std::string refinementID()
	{
		return _refinementID;
	}

	std::string getFilename()
	{
		return _filename;
	}

	std::string metadata()
	{
		return _metadata;
	}

	void setMetadata(std::string metadata)
	{
		_metadata = metadata;
	}
	
	void setPdbPath(std::string path)
	{
		_pdbPath = path;
	}
	
	std::string getPdbPath()
	{
		return _pdbPath;
	}
	
	void setMarked(bool mark)
	{
		if (_dead) return;

		_mark = mark;
		_sele = false;
	}
	
	bool isMarked()
	{
		return _mark;
	}
	
	bool isDead()
	{
		return _dead;
	}
	
	void setDead(bool dead)
	{
		_dead = dead;
		_mark = false;
		_sele = false;
	}
	
	bool isSelected()
	{
		return _sele;
	}
	
	void flipSelected()
	{
		_sele = !_sele;
	}

	void setSelected(bool sele)
	{
		if (!_mark && !_dead)
		{
			_sele = sele;
		}
	}
	
	void recolourVertex(Vertex *v, bool fullDead = false);
	
	void setCrystal(CrystalPtr crystal)
	{
		_crystal = crystal;
	}
	
	CrystalPtr getCrystal()
	{
		return _crystal;
	}
	
	std::vector<vec3> getAtomPositions()
	{
		return _positions;
	}
	
	void setAtomPositions(std::vector<vec3> positions)
	{
		_positions = positions;
	}
	
	void setRWork(double rw)
	{
		_rWork = rw;
	}
	
	void setRFree(double rf)
	{
		_rFree = rf;
	}
	
	double getRWork()
	{
		return _rWork;
	}
	
	double getRFree()
	{
		return _rFree;
	}
private:
	std::string _filename;
	std::string _pdbPath;
	std::string _metadata;
	std::string _refinementID;
	CrystalPtr _crystal;
	std::vector<vec3> _positions;

	bool _mark;
	bool _sele;
	bool _dead;
	
	double _rWork;
	double _rFree;
};

#endif

