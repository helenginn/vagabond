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
	SelectResidue,
	SelectConformer,
} SelectionType;

class Selected2GL : public Atoms2GL
{
public:
	Selected2GL();
	
	void setPicked(AtomPtr atom);
	bool isRefinable();
	AtomGroupPtr refinableSelection();
	void manualRefine();
	void focusOnGroup();
	void splitSelected();
	
	bool isRefining()
	{
		return _refining;
	}
	
	void setMouseRefinement(bool val);
	
	void cancelRefine()
	{
		_refining = false;
	}
	
	AtomGroupPtr getSelected();
	
	void setMouseRay(vec3 ray);
protected:
	virtual int processMolecule(MoleculePtr molecule);
	virtual void bindTextures();
	virtual void findAtoms();

private:
	vec3 averageModelPos();
	AtomPtr _picked;
	SelectionType _sType;
	int _conf;
	vec3 _ray;

	bool _refining;
	bool _mousey;
	std::mutex _switch;
};
