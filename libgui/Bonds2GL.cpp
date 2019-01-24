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

#include "../libsrc/Atom.h"
#include "../libsrc/Molecule.h"
#include "Bonds2GL.h"
#include "Shaders/InkBond_vsh.h"
#include "Shaders/InkBond_fsh.h"

Bonds2GL::Bonds2GL(int average)
{
	_average = average;
	setupAverage();
}

void Bonds2GL::setupAverage()
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

void Bonds2GL::getPositions(AtomPtr minAtom, AtomPtr majAtom, 
                            std::vector<vec3> *min,
                            std::vector<vec3> *maj)
{
	ModelPtr minBond = minAtom->getModel();
	ModelPtr majBond = majAtom->getModel();
	
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

void Bonds2GL::updateModel(int *v, int total, std::vector<vec3> &maj, 
                           std::vector<vec3> &min)
{
	for (int k = 0; k < total; k++)
	{
		if (maj.size() <= k || min.size() <= k)
		{
			continue;
		}

		vec3 majStart = maj[k];
		vec3 minStart = min[k];

		_vertices[*v].pos[0] = majStart.x;
		_vertices[*v].pos[1] = majStart.y;
		_vertices[*v].pos[2] = majStart.z;

		/* Middle two vertices should be the same */
		for (int i = 1; i < 3; i++)
		{
			_vertices[*v + i].pos[0] = (minStart.x + majStart.x) / 2;
			_vertices[*v + i].pos[1] = (minStart.y + majStart.y) / 2;
			_vertices[*v + i].pos[2] = (minStart.z + majStart.z) / 2;
		}

		_vertices[*v+3].pos[0] = minStart.x;
		_vertices[*v+3].pos[1] = minStart.y;
		_vertices[*v+3].pos[2] = minStart.z;

		if (_average)
		{
			for (int i = 0; i < 4; i++)
			{
				_vertices[*v+i].extra[0] = k;
				_vertices[*v+i].extra[1] = 0;

				if (i == 1 || i == 2)
				{
					_vertices[*v + i].extra[1] = 1;
				}
			}

			memcpy(_vertices[*v].normal, _vertices[*v+3].pos,
			       sizeof(GLfloat) * 3);
			memcpy(_vertices[*v+1].normal, _vertices[*v+3].pos,
			       sizeof(GLfloat) * 3);
			memcpy(_vertices[*v+2].normal, _vertices[*v].pos,
			       sizeof(GLfloat) * 3);
			memcpy(_vertices[*v+3].normal, _vertices[*v].pos,
			       sizeof(GLfloat) * 3);

		}

		*v += 4;
	}
}

void Bonds2GL::updateAtoms()
{
	for (int i = 0; i < _pairList.size(); i++)
	{
		AtomWkr wmin = _pairList[i].min;
		AtomWkr wmaj = _pairList[i].maj;
		
		if (wmaj.expired() || wmin.expired())
		{
			continue;
		}

		AtomPtr maj = wmaj.lock();
		AtomPtr min = wmin.lock();

		int total = _pairList[i].size;
		int v = _pairList[i].vNum;

		if (!min->getModel()->isBond()) continue;
		if (!maj->getModel()->isBond()) continue;
		
		std::vector<vec3> majBonds, minBonds;
		getPositions(maj, min, &minBonds, &majBonds);

		updateModel(&v, total, minBonds, majBonds);
	}
}

bool Bonds2GL::addToModel(AtomPtr minor, AtomPtr major, GLuint *count)
{
	std::vector<vec3> minBonds;
	std::vector<vec3> majBonds;
	int start = *count;

	getPositions(minor, major, &minBonds, &majBonds);

	if (!minBonds.size() || !majBonds.size())
	{
		return false;
	}
	
	for (int j = 0; j < majBonds.size(); j++)
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
		setVertexColour(minor, &vertex);
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
			_indices.push_back(*count);
			_indices.push_back(*count + 1);
			_indices.push_back(*count + 2);
			_indices.push_back(*count + 3);
		}
		else if (j == 0)
		{
			/* Suitable for GL_TRIANGLES */
			_indices.push_back(*count + 1);
			_indices.push_back(*count + 0);
			_indices.push_back(*count + 4);

			_indices.push_back(*count + 1);
			_indices.push_back(*count + 4);
			_indices.push_back(*count + 5);

			_indices.push_back(*count + 2);
			_indices.push_back(*count + 3);
			_indices.push_back(*count + 6);

			_indices.push_back(*count + 3);
			_indices.push_back(*count + 6);
			_indices.push_back(*count + 7);
		}

		if (j == 0)
		{
			Atom3D pair;
			pair.maj = major;
			pair.min = minor;
			pair.size = majBonds.size();
			pair.vNum = *count;
			_pairList.push_back(pair);
		}

		*count += 4;
	}

	/* Set up texture, only if using _averages */

	if (!_average)
	{
		return true;
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

	return true;
}

int Bonds2GL::processMolecule(MoleculePtr molecule)
{
	GLuint count = (int)_vertices.size();
	int bonds = 0;

	ExplicitModel::useMutex();

	for (int i = 0; i < molecule->atomCount(); i++)
	{
		AtomPtr minor = molecule->atom(i);

		if (!minor)
		{
			continue;
		}

		if (!minor->getElement())
		{
			continue;
		}

		if (minor->getElectronCount() <= 1)
		{
			continue;
		}

		if (minor->getModel() && !minor->getModel()->isBond())
		{
			continue;
		}

		AtomPtr major = ToBondPtr(minor->getModel())->getMajor();

		if (addToModel(minor, major, &count))
		{
			bonds++;
		}

	}

	_shouldGetBonds = false;

	return bonds;
}

void Bonds2GL::bindTextures()
{
	if (!_average)
	{
		return;
	}

	int num = 1;
	_textures.resize(num);

	glGenTextures(num, &_textures[0]);
	glBindTexture(GL_TEXTURE_2D, _textures[0]);
	checkErrors();
	bindOneTexture(pic_bond);
}

void Bonds2GL::render()
{
	if (!_enabled)
	{
		return;
	}

	Vagabond2GL::render();

	if (_average)
	{
		reorderIndices();
	}

	GLObject::render();


}
