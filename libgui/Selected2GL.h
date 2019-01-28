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

#include "Atoms2GL.h"

typedef enum
{
	SelectAtom,
	SelectSidechain,
	SelectConformer,
	SelectMonomer,
	SelectMonConf,
} SelectionType;

class Selected2GL : public Atoms2GL
{
public:
	Selected2GL();
	
	void setPicked(AtomPtr atom, bool preserveType = false);
	AtomGroupPtr refinableSelection();
	void manualRefine();
	void focusOnGroup();
	void splitSelected();
	void deleteSelected();
	void toggleKicks();
	void advanceMonomer(int dir);
	void addPicked(AtomPtr atom, bool preserveType);

	void setAdding(bool val)
	{
		_adding = val;
	}
	
	bool isRefining()
	{
		return _refining;
	}
	
	void setMouseRefinement(bool val);
	
	void cancelRefine();
	
	AtomPtr getPicked()
	{
		if (_picked.size())
		{
			return _picked.back();
		}

		return AtomPtr();
	}
	
	AtomGroupPtr getSelected();
	
	void setMouseRay(vec3 ray);
protected:
	virtual int processMolecule(MoleculePtr molecule);
	virtual void bindTextures();
	virtual void findAtoms();

private:
	vec3 averageModelPos();
	std::vector<AtomPtr> _picked;
	SelectionType _sType;
	int _conf;
	vec3 _ray;

	bool _adding;
	bool _refining;
	bool _kicking;
	bool _mousey;
	std::mutex _switch;
};
