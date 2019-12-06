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
#include "Blob.h"
#include "../libsrc/Molecule.h"
#include "Shaders/Blob_vsh.h"
#include "Shaders/Blob_fsh.h"

Atoms2GL::Atoms2GL()
{
	_renderType = GL_POINTS;
	_vertShader = &Blob_vsh;
	_fragShader = &Blob_fsh;
	_usesFocalDepth = true;
}

void Atoms2GL::addAtom(AtomPtr atom)
{
	std::vector<vec3> atomPos;
	getPositions(atom, AtomPtr(), &atomPos, NULL);
	vec3 pos = atomPos[0];

	Vertex vertex;
	vertex.pos[0] = pos.x;
	vertex.pos[1] = pos.y;
	vertex.pos[2] = pos.z;
	setVertexColour(atom, &vertex);
	GLuint count = (int)_vertices.size();
	_vertices.push_back(vertex);
	_indices.push_back(count);
	Atom3D pair;
	pair.min = atom;
	pair.size = 0;
	pair.vNum = 0;
	_pairList.push_back(pair);
}

int Atoms2GL::processMolecule(MoleculePtr molecule)
{
	GLuint count = (int)_vertices.size();

	ExplicitModel::useMutex();

	for (int i = 0; i < molecule->atomCount(); i++)
	{
		AtomPtr atom = molecule->atom(i);
		
		if (!atom)
		{
			continue;
		}

		if (!atom->getElement())
		{
			continue;
		}

		if (atom->getElectronCount() <= 1)
		{
			continue;
		}

		/* We only want to draw non-explicit models in this instance */
		if (atom->getModel() && atom->getModel()->hasExplicitPositions())
		{
			continue;
		}

		addAtom(atom);
		count++;
	}
	
	return molecule->atomCount();
	_shouldGetBonds = false;
}

bool Atoms2GL::getPositions(AtomPtr minor, AtomPtr major, 
                            std::vector<vec3> *min,
                            std::vector<vec3> *maj)
{
	vec3 ave = minor->getAbsolutePosition();
	*min = std::vector<vec3>(1, ave);
	
	return true;
}

void Atoms2GL::bindTextures()
{
	Vagabond2GL::bindTextures();
	bindOneTexture(pic_atom_blob);
}


void Atoms2GL::updateAtoms()
{

}

void Atoms2GL::render()
{
	if (!_enabled)
	{
		return;
	}
	
	Vagabond2GL::render();
	GLObject::render();
}

