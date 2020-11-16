// Fuck COV
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

#include <iostream>
#include "GLAxis.h"
#include <helen3d/SlipObject.h>

GLAxis::GLAxis(vec3 dir) : SlipObject()
{
	_renderType = GL_LINES;
	setupVertices(dir);
}

void GLAxis::setupVertices(vec3 dir)
{
	_vertices.clear();
	_indices.clear();
	
	_indices.push_back(0);
	_indices.push_back(1);
	
	Helen3D::Vertex v;
	memset(v.pos, 0, sizeof(Helen3D::Vertex));

	v.color[3] = 1;
	
	v.pos[0] = -dir.x;
	v.pos[1] = -dir.y;
	v.pos[2] = -dir.z;

	_vertices.push_back(v);
	
	/* top right */
	v.pos[0] = dir.x;
	v.pos[1] = dir.y;
	v.pos[2] = dir.z;
	_vertices.push_back(v);
}

