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
#include "../../libsrc/Model.h"
#include "../../libsrc/Atom.h"
#include "Shaders/InkBond_vsh.h"
#include "Shaders/InkBond_fsh.h"
#include "Shaders/SlimBond_fsh.h"
#include "Shaders/SlimBond_vsh.h"
#include "Dotted.h"

using namespace Helen3D;

Connect2GL::Connect2GL() : Vagabond2GL()
{
	_renderType = GL_TRIANGLES;
	_vString = InkBond_vsh();
	_fString = InkBond_fsh();
	setNeedsExtra(true);
	_shouldGetBonds = false;
	_usesFocalDepth = true;
}

bool Connect2GL::getPositions(AtomPtr minor, AtomPtr major, 
                              std::vector<vec3> *min,
                              std::vector<vec3> *maj)
{
	vec3 ave = minor->getAbsolutePosition();
	*min = std::vector<vec3>(1, ave);
	
	ave = major->getAbsolutePosition();
	*maj = std::vector<vec3>(1, ave);
	
	return true;
}

void Connect2GL::bindTextures()
{
	bindOneTexture(pic_dotted_line);
}

void Connect2GL::addAtoms(AtomPtr major, AtomPtr minor)
{
	std::vector<vec3> majPos, minPos;
	getPositions(minor, major, &minPos, &majPos);

	vec3 maj = majPos[0];
	vec3 min = minPos[0];

	/*
	vec3 normal = vec3_subtract_vec3(min, maj);
	vec3_set_length(&normal, 1);
	GLfloat glNorm[3];
	glNorm[0] = normal.x;
	glNorm[1] = normal.y;
	glNorm[2] = normal.z;
	*/
	
	int n = _vertices.size();
	Vertex vertex;
	memset(&vertex, '\0', sizeof(Vertex));
	setVertexColour(minor, &vertex);

	/* major */
	{
		vertex.pos[0] = maj.x;
		vertex.pos[1] = maj.y;
		vertex.pos[2] = maj.z;
		vertex.tex[0] = 0;
		vertex.tex[1] = 0;
		vertex.extra[0] = 0;
		_vertices.push_back(vertex);
		vertex.tex[0] = 0.0;
		vertex.tex[1] = 1.0;
		vertex.extra[0] = 1;
		_vertices.push_back(vertex);
	}

	/* minor */
	{
		vertex.pos[0] = min.x;
		vertex.pos[1] = min.y;
		vertex.pos[2] = min.z;
		vertex.tex[0] = 1.0;
		vertex.tex[1] = 0.0;
		vertex.extra[0] = 1;
		_vertices.push_back(vertex);
		vertex.tex[0] = 1.0;
		vertex.tex[1] = 1.0;
		vertex.extra[0] = 0;
		_vertices.push_back(vertex);
	}
	
	memcpy(_vertices[n].normal, _vertices[n+2].pos, sizeof(GLfloat) * 3);
	memcpy(_vertices[n+1].normal, _vertices[n+3].pos, sizeof(GLfloat) * 3);
	memcpy(_vertices[n+2].normal, _vertices[n+0].pos, sizeof(GLfloat) * 3);
	memcpy(_vertices[n+3].normal, _vertices[n+1].pos, sizeof(GLfloat) * 3);

	_indices.push_back(n + 0);
	_indices.push_back(n + 1);
	_indices.push_back(n + 2);

	_indices.push_back(n + 1);
	_indices.push_back(n + 2);
	_indices.push_back(n + 3);

	Atom3D pair;
	pair.maj = major;
	pair.min = minor;
	pair.size = 1;
	pair.vNum = n;
	_pairList.push_back(pair);
}

void Connect2GL::updateAtoms()
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

		int v = _pairList[i].vNum;
		std::vector<vec3> majPos, minPos;
		getPositions(min, maj, &minPos, &majPos);

		for (int j = v; j <= v + 1; j++)
		{
			_vertices[j].pos[0] = majPos[0].x;
			_vertices[j].pos[1] = majPos[0].y;
			_vertices[j].pos[2] = majPos[0].z;
		}

		for (int j = v + 2; j <= v + 3; j++)
		{
			_vertices[j].pos[0] = minPos[0].x;
			_vertices[j].pos[1] = minPos[0].y;
			_vertices[j].pos[2] = minPos[0].z;
		}
	}
}

void Connect2GL::clear()
{
	_pairList.clear();
	_vertices.clear();
	_indices.clear();
}

void Connect2GL::render(SlipGL *sender)
{
	if (!_enabled)
	{
		return;
	}
	
	Vagabond2GL::render(sender);
	reorderIndices();
	SlipObject::render(sender);
}
