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

#include "Vagabond2GL.h"

class Bonds2GL : public Vagabond2GL
{
public:
	Bonds2GL(int average = false);

	virtual void render();

protected:
	virtual void updateAtoms();
	/** returns whether bonds could be fished */
	virtual bool getPositions(AtomPtr minAtom, AtomPtr majAtom, 
	                          std::vector<vec3> *min,
	                          std::vector<vec3> *maj);

	bool addToModel(AtomPtr minor, AtomPtr major, GLuint *count);

	virtual int processMolecule(MoleculePtr molecule);
	virtual void bindTextures();
private:

	virtual bool acceptablePositions(AtomPtr minAtom);
	void setupAverage();
	void updateModel(int *v, int total, std::vector<vec3> &maj, 
	                 std::vector<vec3> &min, AtomPtr atm);

	int _average;
};
