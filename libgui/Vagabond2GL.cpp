//
//  Vagabond2GL.cpp
//  VagabondViewer
//
//  Created by Helen Ginn on 03/11/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#include "Vagabond2GL.h"
#include "Options.h"
#include "Crystal.h"
#include "Bond.h"
#include "Atom.h"
#include "Element.h"

void Vagabond2GL::updateAtoms()
{
	for (AtomMap::iterator it = _atomMap.begin(); it != _atomMap.end(); it++)
	{
		AtomPtr atom = it->first;
		std::pair<int, int> pair = it->second;
		int conformer = pair.first;
		int vertex = pair.second;

		if (!atom->getModel()->isBond()) continue;

		std::vector<vec3> majBonds, minBonds;
		getPositions(atom, &minBonds, &majBonds);

		if (majBonds.size() >= conformer || minBonds.size() >= conformer)
		{
			continue;
		}

		vec3 majStart = majBonds[conformer];
		vec3 minStart = minBonds[conformer];

		_vertices[vertex].pos[0] = majStart.x;
		_vertices[vertex].pos[1] = majStart.y;
		_vertices[vertex].pos[2] = majStart.z;
		vertex++;
		_vertices[vertex].pos[0] = minStart.x;
		_vertices[vertex].pos[1] = minStart.y;
		_vertices[vertex].pos[2] = minStart.z;

	}
}

void Vagabond2GL::getPositions(AtomPtr atom, std::vector<vec3> *min,
							   std::vector<vec3> *maj)
{
	ModelPtr minBond = atom->getModel();
	ModelPtr majBond = (ToBondPtr(minBond))->getMajor()->getModel();

	if (majBond->isBond())
	{
		*maj = ToBondPtr(majBond)->fishPositions();
	}

	*min = ToBondPtr(minBond)->fishPositions();
}

bool Vagabond2GL::shouldGetBonds()
{
	if (!_moleculeMap.size())
	{
		return true;
	}

	_renders = 0;

	OptionsPtr globalOptions = Options::getRuntimeOptions();

	for (int i = 0; i < globalOptions->crystalCount(); i++)
	{
		CrystalPtr crystal = globalOptions->getCrystal(i);

		for (int j = 0; j < crystal->moleculeCount(); j++)
		{
			MoleculePtr molecule = crystal->molecule(j);

			if (!_moleculeMap.count(molecule))
			{
				return true;
			}

			int expected = _moleculeMap[molecule];

			int existing = 0;

			for (int k = 0; k < molecule->atomCount(); k++)
			{
				if (molecule->atom(k)->getModel()->isBond())
				{
					existing++;
				}
			}

			if (expected != existing)
			{
				return true;
			}
		}
	}

	return false;
}

int Vagabond2GL::processMolecule(MoleculePtr molecule)
{
	GLuint count = (int)_vertices.size();
	int bonds = 0;

	for (int i = 0; i < molecule->atomCount(); i++)
	{
		AtomPtr atom = molecule->atom(i);

		if (atom->getModel()->isBond())
		{
			ToBondPtr(atom->getModel())->useMutex();
		}
	}

	for (int i = 0; i < molecule->atomCount(); i++)
	{
		AtomPtr atom = molecule->atom(i);

		if (atom->getElement()->electronCount() <= 1)
		{
			continue;
		}

		if (atom->getModel()->isBond())
		{
			std::vector<vec3> majBonds, minBonds;
			getPositions(atom, &minBonds, &majBonds);

			for (int j = 0; j < majBonds.size(); j += 1)
			{
				vec3 majStart = majBonds[j];
				if (majBonds.size() <= j || minBonds.size() <= j)
				{
					break;
				}

				vec3 minStart = minBonds[j];
				Vertex vertex;
				vertex.pos[0] = majStart.x;
				vertex.pos[1] = majStart.y;
				vertex.pos[2] = majStart.z;
				vertex.color[0] = 0;
				vertex.color[1] = 0;
				vertex.color[2] = 0;
				vertex.color[3] = 1;
				_vertices.push_back(vertex);
				vertex.pos[0] = minStart.x;
				vertex.pos[1] = minStart.y;
				vertex.pos[2] = minStart.z;
				_vertices.push_back(vertex);
				_indices.push_back(count);
				_indices.push_back(count + 1);

				_atomMap[atom] = std::make_pair(j, count);
				count += 2;
			}

			if (minBonds.size() && majBonds.size())
			{
				bonds++;
			}
		}
	}

	return bonds;
}

void Vagabond2GL::findAtoms()
{
	OptionsPtr globalOptions = Options::getRuntimeOptions();

	_vertices.clear();
	_indices.clear();
	_atomMap.clear();
	_moleculeMap.clear();

	if (globalOptions)
	{
		for (int i = 0; i < globalOptions->crystalCount(); i++)
		{
			CrystalPtr crystal = globalOptions->getCrystal(i);

			for (int j = 0; j < crystal->moleculeCount(); j++)
			{
				MoleculePtr molecule = crystal->molecule(j);

				int expected = processMolecule(molecule);

				_moleculeMap[molecule] = expected;
			}
		}
	}
}

void Vagabond2GL::render()
{
	if (shouldGetBonds())
	{
		findAtoms();
		rebindProgram();
	}
	else
	{
		updateAtoms();
	}

	GLObject::render();
}
