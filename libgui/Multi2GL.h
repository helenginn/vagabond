// Vagabond
// Copyright (C) 2020 Helen Ginn
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

#ifndef __vagabond__multi2gl__
#define __vagabond__multi2gl__

#include "Atoms2GL.h"
#include "../../libsrc/shared_ptrs.h"

class Connect2GL;

class Multi2GL : public Atoms2GL
{
public:
	Multi2GL();

	virtual void render(SlipGL *sender);
	
	Connect2GLPtr getConnected2GL()
	{
		return _connected;
	}
protected:
	virtual int processMolecule(MoleculePtr molecule);
	virtual void updateAtoms();
	virtual void bindTextures();
	virtual bool getPositions(AtomPtr major, AtomPtr minor,
	                          std::vector<vec3> *min,
	                          std::vector<vec3> *maj);

	void addAtom(AtomPtr atom);
	void addConnections(SpongePtr sponge);
private:
	Connect2GLPtr _connected;

};

#endif
