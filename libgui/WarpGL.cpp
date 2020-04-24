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

#include "WarpGL.h"
#include "Arrow.h"
#include "../libsrc/Crystal.h"
#include "../libsrc/Options.h"
#include "../libsrc/SpaceWarp.h"
#include "Shaders/PointyTriangle_vsh.h"
#include "Shaders/InkBond_fsh.h"

WarpGL::WarpGL() : GLObject()
{
	_renderType = GL_TRIANGLES;
	_vertShader = PointyTriangle_vsh();
	_fragShader = InkBond_fsh();
	_extra = true;
	_usesFocalDepth = true;
}

void WarpGL::initialise()
{
	CrystalPtr crystal = Options::getActiveCrystal();
	
	if (!crystal)
	{
		return;
	}
	
	_warp = crystal->getWarp();
	
	if (!_warp)
	{
		return;
	}
	
	_fft = _warp->getFFT();
	
	if (!_fft)
	{
		return;
	}
	
	generateVertices();
}

inline void closestToFocus(vec3 *p, vec3 f)
{
	while (p->x - f.x < -0.5) p->x++;
	while (p->x - f.x >= 0.6) p->x--;

	while (p->y - f.y < -0.5) p->y++;
	while (p->y - f.y >= 0.6) p->y--;

	while (p->z - f.z < -0.5) p->z++;
	while (p->z - f.z >= 0.6) p->z--;

}

void WarpGL::wrapVerticesRound()
{
	if (!_fft)
	{
		return;
	}
	
	mat4x4 inv = mat4x4_inverse(modelMat);
	mat4x4 trans = mat4x4_transpose(inv);
	
	mat3x3 toReal = _fft->toReal();
	mat3x3 toRecip = _fft->toRecip();

	vec3 f = make_vec3(0, 0, -13);
	f = mat4x4_mult_vec(inv, f);
	mat3x3_mult_vec(toRecip, &f);

	for (int i = 0; i < _vertices.size(); i++)
	{
		Vertex *v = &_vertices[i];
		vec3 p = make_vec3(v->pos[0], v->pos[1], v->pos[2]);
		mat3x3_mult_vec(toRecip, &p);
		closestToFocus(&p, f);
		mat3x3_mult_vec(toReal, &p);
		vec3ToVertex(*v, p);
	}
}

void WarpGL::generateVertices()
{
	/* position = position of vector,
	 * normal = direction vector,
	 * tex = (-1, -1), (-1, 1), (1, -1) or (1, 1) */

	long count = 0;
	mat3x3 toReal = _fft->toReal();

	for (long i = 0; i < _fft->nn(); i++)
	{
		vec3 offset = _warp->getWarp(i);
		
		if (vec3_length(offset) < 1e-6)
		{
			/* not going to display a zero */
			continue;
		}

		Vertex v;
		memset(&v, '\0', sizeof(Vertex));

		vec3 pos = _fft->fracFromElement(i);
		mat3x3_mult_vec(toReal, &pos);
		vec3ToVertex(v, pos);
		vec3ToNormal(v, offset);
		v.color[0] = 0.5;
		v.color[1] = 0.5;
		v.color[2] = 1.0;
		v.color[3] = 1.0;
		
		_indices.push_back(count + 0);
		_indices.push_back(count + 1);
		_indices.push_back(count + 2);
		_indices.push_back(count + 2);
		_indices.push_back(count + 1);
		_indices.push_back(count + 3);

		for (int m = -1; m <= 1; m+=2)
		{
			for (int n = -1; n <= 1; n+=2)
			{
				v.tex[0] = (m < 0 ? 0 : m);
				v.tex[1] = (n < 0 ? 0 : n);
				v.extra[0] = m;
				v.extra[1] = n;
				
				_vertices.push_back(v);
				count++;
			}
		}
	}
}

void WarpGL::render()
{
	if (_vertices.size() == 0)
	{
		initialise();
	}
	
	wrapVerticesRound();
	
	rebindProgram();
	GLObject::render();
}

void WarpGL::refresh()
{
	_vertices.clear();
	_indices.clear();
}

WarpGL::~WarpGL()
{

}

void WarpGL::bindTextures()
{
	int num = 1;
	_textures.resize(num);

	glGenTextures(num, &_textures[0]);
	glBindTexture(GL_TEXTURE_2D, _textures[0]);
	checkErrors();

	bindOneTexture(pic_arrow);
}

