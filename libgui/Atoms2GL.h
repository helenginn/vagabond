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

#ifndef __vagabond__atoms2gl__
#define __vagabond__atoms2gl__

#include "Vagabond2GL.h"

class Atoms2GL : public Vagabond2GL
{
public:
	Atoms2GL();

	virtual void render(SlipGL *sender);
protected:
	virtual int processMolecule(MoleculePtr molecule);
	virtual void updateAtoms();
	virtual void bindTextures();
	virtual bool getPositions(AtomPtr major, AtomPtr minor,
	                          std::vector<vec3> *min,
	                          std::vector<vec3> *maj);

	void addAtom(AtomPtr atom);
private:

};

#endif
