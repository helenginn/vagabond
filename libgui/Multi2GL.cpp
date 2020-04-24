// Vagabond
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

#include "Connect2GL.h"
#include "Multi2GL.h"
#include "Cross.h"
#include "Shaders/Blob_vsh.h"
#include "Shaders/Blob_fsh.h"
#include "../../libsrc/Molecule.h"
#include "../../libsrc/Model.h"
#include "../../libsrc/Atom.h"
#include "../../libsrc/Sponge.h"

Multi2GL::Multi2GL()
{
	_renderType = GL_POINTS;
	_vertShader = Blob_vsh();
	_fragShader = Blob_fsh();
	_usesFocalDepth = true;
	
	_connected = Connect2GLPtr(new Connect2GL());
}

void Multi2GL::addAtom(AtomPtr atom)
{
	if (!atom->getModel()->hasExplicitPositions())
	{
		return;
	}

	std::vector<vec3> atomPos;
	getPositions(atom, AtomPtr(), &atomPos, NULL);
	GLuint count = (int)_vertices.size();
	
	if (atomPos.size() == 0)
	{
		return;
	}

	Atom3D pair;
	pair.min = atom;
	pair.size = atomPos.size();
	pair.vNum = count;
	_pairList.push_back(pair);

	for (int i = 0; i < atomPos.size(); i++)
	{
		vec3 pos = atomPos[i];
		Vertex vertex;
		vertex.pos[0] = pos.x;
		vertex.pos[1] = pos.y;
		vertex.pos[2] = pos.z;
		vertex.normal[0] = 0.25;
		setVertexColour(atom, &vertex);
		_vertices.push_back(vertex);
		_indices.push_back(count);
		count++;
	}
	
	if (atom->getModel()->isSponge())
	{
		addConnections(ToSpongePtr(atom->getModel()));
	}
}

int Multi2GL::processMolecule(MoleculePtr molecule)
{
	GLuint count = (int)_vertices.size();

	ExplicitModel::useMutex();
	
	_connected->clear();

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

		/* We only want to draw explicit models in this instance */
		if (atom->getModel() && !atom->getModel()->hasExplicitPositions())
		{
			continue;
		}

		if (atom->getModel()->isAnchor() || atom->getModel()->isBond())
		{
			continue;
		}

		addAtom(atom);
		count++;
	}
	
	_shouldGetBonds = false;
	return molecule->atomCount();
}

void Multi2GL::addConnections(SpongePtr sponge)
{
	AtomPtr a = sponge->getAtom(); 
	for (int i = 0; i < sponge->connectedCount(); i++)
	{
		AtomPtr c = sponge->getConnectedAtom(i);
		_connected->addAtoms(a, c);
	}
}

bool Multi2GL::getPositions(AtomPtr minor, AtomPtr major, 
                            std::vector<vec3> *min,
                            std::vector<vec3> *maj)
{
	ExplicitModelPtr minModel = minor->getExplicitModel();
	vec3 minAve;
	*min = minModel->fishPositions(&minAve);
	
	if (min->size() != _lastEnsembleCount && min->size() >= 0)
	{
		_lastEnsembleCount = min->size();
		_shouldGetBonds = true;
	}
	
	return true;
}

void Multi2GL::bindTextures()
{
	Vagabond2GL::bindTextures();
	bindOneTexture(pic_atom_cross);
}

void Multi2GL::updateAtoms()
{
	for (int i = 0; i < _pairList.size(); i++)
	{
		Atom3D entry = _pairList[i];
		
		std::vector<vec3> atomPos;
		getPositions(entry.min.lock(), AtomPtr(), &atomPos, NULL);

		int count = 0;
		for (int j = 0; j < atomPos.size(); j++)
		{
			vec3 pos = atomPos[j];
			Vertex *v = &_vertices[entry.vNum + count];
			v->pos[0] = pos.x;
			v->pos[1] = pos.y;
			v->pos[2] = pos.z;
			count++;
		}
	}
}

void Multi2GL::render()
{
	Atoms2GL::render();
}
