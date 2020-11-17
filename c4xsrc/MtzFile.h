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

#include <shared_ptrs.h>
#include <vec3.h>
#include <string>
#include <mat3x3.h>
#include <Frameworks.h>

struct _Vertex;
typedef _Vertex Vertex;

class QuickAtoms;

class MtzFile
{
public:
	MtzFile(std::string filename);
	
	void setRefinementID(std::string id)
	{
		_refinementID = id;
	}
	
	void setPanddaName(std::string file)
	{
		_panddaName = file;
	}
	
	void setCifPath(std::string file)
	{
		_cifPath = file;
	}
	
	std::string getCifPath()
	{
		return _cifPath;
	}
	
	void setLigPath(std::string file)
	{
		_ligPath = file;
	}
	
	std::string getLigPath()
	{
		return _ligPath;
	}
	
	std::string refinementID()
	{
		return _refinementID;
	}

	std::string getFilename()
	{
		return _filename;
	}

	std::string getPanddaName()
	{
		return _panddaName;
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
	
	void recolourVertex(Helen3D::Vertex *v, bool fullDead = false);
	
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
	
	QuickAtoms *getQuickAtoms()
	{
		return _quickAtoms;
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
	
	void setColour(double r, double g, double b, double a);
private:
	std::string _filename;
	std::string _pdbPath;
	std::string _cifPath;
	std::string _ligPath;
	std::string _metadata;
	std::string _refinementID;
	std::string _panddaName;
	CrystalPtr _crystal;
	std::vector<vec3> _positions;
	QuickAtoms *_quickAtoms;

	bool _mark;
	bool _sele;
	bool _dead;
	
	double _rWork;
	double _rFree;
	
	double _red;
	double _green;
	double _blue;
	double _alpha;
};

#endif

