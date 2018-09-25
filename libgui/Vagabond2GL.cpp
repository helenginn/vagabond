//
//  Vagabond2GL.cpp
//  VagabondViewer
//
//  Created by Helen Ginn on 03/11/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#include "Vagabond2GL.h"
#include "../libsrc/Options.h"
#include "../libsrc/Crystal.h"
#include "../libsrc/Bond.h"
#include "../libsrc/Atom.h"
#include "../libsrc/Element.h"
#include "Shaders/InkBond_vsh.h"
#include "Shaders/InkBond_fsh.h"

void Vagabond2GL::setupAverage()
{
	if (!_average)
	{
		return;
	}

	_renderType = GL_TRIANGLES;
	_vertShader = &InkBond_vsh;
	_fragShader = &InkBond_fsh;
	_extra = true;
}

void Vagabond2GL::updateAtoms()
{
	for (AtomMap::iterator it = _atomMap.begin(); it != _atomMap.end(); it++)
	{
		AtomPtr atom = it->first;
		std::pair<int, int> pair = it->second;
		int total = pair.first;
		int v = pair.second;
		//		std::cout << "c/v: " << conformer << " " << v << std::endl;

		if (!atom->getModel()->isBond()) continue;

		std::vector<vec3> majBonds, minBonds;
		getPositions(atom, &minBonds, &majBonds);

		for (int k = 0; k < total; k++)
		{
			if (majBonds.size() <= k || minBonds.size() <= k)
			{
				continue;
			}

			vec3 majStart = majBonds[k];
			vec3 minStart = minBonds[k];

			_vertices[v].pos[0] = majStart.x;
			_vertices[v].pos[1] = majStart.y;
			_vertices[v].pos[2] = majStart.z;

			/* Middle two vertices should be the same */
			for (int i = 1; i < 3; i++)
			{
				_vertices[v + i].pos[0] = (minStart.x + majStart.x) / 2;
				_vertices[v + i].pos[1] = (minStart.y + majStart.y) / 2;
				_vertices[v + i].pos[2] = (minStart.z + majStart.z) / 2;
			}

			_vertices[v+3].pos[0] = minStart.x;
			_vertices[v+3].pos[1] = minStart.y;
			_vertices[v+3].pos[2] = minStart.z;

			if (_average)
			{
				for (int i = 0; i < 4; i++)
				{
					_vertices[v+i].extra[0] = k;
					_vertices[v+i].extra[1] = 0;

					if (i == 1 || i == 2)
					{
						_vertices[v + i].extra[1] = 1;
					}
				}

				memcpy(_vertices[v].normal, _vertices[v+3].pos,
				       sizeof(GLfloat) * 3);
				memcpy(_vertices[v+1].normal, _vertices[v+3].pos,
				       sizeof(GLfloat) * 3);
				memcpy(_vertices[v+2].normal, _vertices[v].pos,
				       sizeof(GLfloat) * 3);
				memcpy(_vertices[v+3].normal, _vertices[v].pos,
				       sizeof(GLfloat) * 3);
				
				for (int i = 0; false && i < 4; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						std::cout << _vertices[v + i].pos[j] << " ";
					}

					std::cout << " - ";

					for (int j = 0; j < 3; j++)
					{
						std::cout << _vertices[v + i].normal[j] << " ";
					}

					std::cout << "(" << _vertices[v + i].extra[0] << ")";

					std::cout << std::endl;
				}
			}
			
			v += 4;
		}
	}
}

void Vagabond2GL::getPositions(AtomPtr atom, std::vector<vec3> *min,
                               std::vector<vec3> *maj)
{
	ModelPtr minBond = atom->getModel();
	ModelPtr majBond = (ToBondPtr(minBond))->getMajor()->getModel();
	
	vec3 minAve, majAve;

	*maj = ToBondPtr(majBond)->fishPositions(&majAve);
	*min = ToBondPtr(minBond)->fishPositions(&minAve);
	
	if (_average)
	{
		*min = std::vector<vec3>(2, minAve);
		*maj = std::vector<vec3>(2, majAve);
	}
	
	if (maj->size() != _lastEnsembleCount)
	{
		_lastEnsembleCount = maj->size();
		_shouldGetBonds = true;
	}
}

bool Vagabond2GL::shouldGetBonds()
{
	if (!_moleculeMap.size() || _shouldGetBonds)
	{
		return true;
	}

	_renders = 0;
	
	/* Have there been any changes? */

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
				AtomPtr thisAtom = molecule->atom(k);
				
				if (!thisAtom)
				{
					continue;
				}

				if (!thisAtom->getModel() ||
				    !thisAtom->getModel()->isBond())
				{
					continue;
				}
				
				if (!thisAtom->getElement() || 
				    thisAtom->getElement()->electronCount() <= 1)
				{
					continue;
				}

				existing++;
			}

			if (expected != existing)
			{
				return true;
			}
		}
	}

	return false;
}

void Vagabond2GL::setVertexColour(AtomPtr atom, Vertex *vertex)
{
	vertex->color[0] = 100. / 255.;
	vertex->color[1] = 100. / 255.;
	vertex->color[2] = 100. / 255.;
	vertex->color[3] = 1.0;

	if (atom->getElement()->getSymbol() == "O")
	{
		vertex->color[0] = 1.0;
		vertex->color[1] = 0.0;
		vertex->color[2] = 0.0;
	}

	if (atom->getElement()->getSymbol() == "N")
	{
		vertex->color[0] = 104. / 255.;
		vertex->color[1] = 139. / 255.;
		vertex->color[2] = 255. / 255.;
	}

	if (atom->getElement()->getSymbol() == "S")
	{
		vertex->color[0] = 255. / 255.;
		vertex->color[1] = 255. / 255.;
		vertex->color[2] = 0. / 255.;
	}
}

int Vagabond2GL::processMolecule(MoleculePtr molecule)
{
	GLuint count = (int)_vertices.size();
	int bonds = 0;

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

		if (atom->getElement()->electronCount() <= 1)
		{
			continue;
		}


		if (atom->getModel() && !atom->getModel()->isBond())
		{
			continue;
		}

		if (atom->getModel()->isBond())
		{
			AtomPtr major = ToBondPtr(atom->getModel())->getMajor();

			std::vector<vec3> majBonds, minBonds;
			getPositions(atom, &minBonds, &majBonds);
			
			int start = count;

			for (int j = 0; j < majBonds.size(); j += 1)
			{
				vec3 majStart = majBonds[j];
				if (majBonds.size() <= j || minBonds.size() <= j)
				{
					break;
				}


				vec3 minStart = minBonds[j];
				vec3 normal = vec3_subtract_vec3(minStart, majStart);
				vec3_set_length(&normal, 1);
				GLfloat glNorm[3];
				glNorm[0] = normal.x;
				glNorm[1] = normal.y;
				glNorm[2] = normal.z;

				Vertex vertex;
				vertex.pos[0] = majStart.x;
				vertex.pos[1] = majStart.y;
				vertex.pos[2] = majStart.z;
				memcpy(vertex.normal, &glNorm, 3 * sizeof(GLfloat));
				setVertexColour(major, &vertex);
				_vertices.push_back(vertex);

				/* Mid point */
				vertex.pos[0] = (minStart.x + majStart.x) / 2;
				vertex.pos[1] = (minStart.y + majStart.y) / 2;
				vertex.pos[2] = (minStart.z + majStart.z) / 2;
				memcpy(vertex.normal, &glNorm, 3 * sizeof(GLfloat));
				_vertices.push_back(vertex);
				setVertexColour(atom, &vertex);
				_vertices.push_back(vertex);

				/* Minor atom */
				vertex.pos[0] = minStart.x;
				vertex.pos[1] = minStart.y;
				vertex.pos[2] = minStart.z;
				memcpy(vertex.normal, &glNorm, 3 * sizeof(GLfloat));
				_vertices.push_back(vertex);
				
				if (_renderType == GL_LINES)
				{
					/* Suitable for GL_LINES */
					_indices.push_back(count);
					_indices.push_back(count + 1);
					_indices.push_back(count + 2);
					_indices.push_back(count + 3);
				}
				else if (j == 0)
				{
					/* Suitable for GL_TRIANGLES */
					_indices.push_back(count + 1);
					_indices.push_back(count + 0);
					_indices.push_back(count + 4);

					_indices.push_back(count + 1);
					_indices.push_back(count + 4);
					_indices.push_back(count + 5);

					_indices.push_back(count + 2);
					_indices.push_back(count + 3);
					_indices.push_back(count + 6);

					_indices.push_back(count + 3);
					_indices.push_back(count + 6);
					_indices.push_back(count + 7);
				}

				if (j == 0)
				{
					_atomMap[atom] = std::make_pair(majBonds.size(), count);
				}

				count += 4;
			}
			
			if (minBonds.size() && majBonds.size())
			{
				bonds++;
			}
			
			/* Set up texture, only if using _averages */

			if (!_average)
			{
				continue;
			}
			
			_vertices[start].tex[0] = 0;
			_vertices[start].tex[1] = 0;

			_vertices[start+1].tex[0] = 0.5;
			_vertices[start+1].tex[1] = 0.0;

			_vertices[start+2].tex[0] = 0.5;
			_vertices[start+2].tex[1] = 0.0;

			_vertices[start+3].tex[0] = 1.0;
			_vertices[start+3].tex[1] = 0.0;

			_vertices[start+4].tex[0] = 0;
			_vertices[start+4].tex[1] = 1.0;

			_vertices[start+5].tex[0] = 0.5;
			_vertices[start+5].tex[1] = 1.0;

			_vertices[start+6].tex[0] = 0.5;
			_vertices[start+6].tex[1] = 1.0;

			_vertices[start+7].tex[0] = 1.0;
			_vertices[start+7].tex[1] = 1.0;
		}
	}
	
	_shouldGetBonds = false;
	return bonds;
}

void Vagabond2GL::bindTextures()
{
	int num = 1;
	_textures.resize(num);

	glGenTextures(num, &_textures[0]);
	glBindTexture(GL_TEXTURE_2D, _textures[0]);
	checkErrors();
	bindOneTexture(pic_bond);
}

void Vagabond2GL::findAtoms()
{
	OptionsPtr globalOptions = Options::getRuntimeOptions();

	clearVertices();
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
	if (!_enabled)
	{
		return;
	}

	if (shouldGetBonds())
	{
		findAtoms();
	}
	else
	{
		updateAtoms();
	}

	rebindProgram();
	
	if (_average)
	{
		reorderIndices();
	}
	GLObject::render();
}
