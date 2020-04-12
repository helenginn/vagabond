// Fuck COV
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

#include <iostream>
#include "GLPoint.h"
#include "Averager.h"
#include "KeeperGL.h"
#include "MtzFFT.h"
#include "MtzFile.h"
#include "shaders/Blob_vsh.h"
#include "shaders/Blob_fsh.h"

GLPoint::GLPoint() : GLObject()
{
	_renderType = GL_POINTS;
	_a = 0;
	_b = 1;
	_c = 2;
}

void GLPoint::initialisePrograms()
{
	GLObject::initialisePrograms(&Blob_vsh, &Blob_fsh);
}

void GLPoint::addPoint(vec3 point)
{
	_indices.push_back(_vertices.size());

	Vertex v;
	memset(v.pos, 0, sizeof(Vertex));

	v.color[3] = 1;
	
	v.pos[0] = point.x;
	v.pos[1] = point.y;
	v.pos[2] = point.z;

	_vertices.push_back(v);
}

void GLPoint::selectInWindow(float x1, float y1, float x2, float y2,
                             int add)
{
	mat4x4 model = _keeper->getModel();
	
	for (size_t i = 0; i < _ave->mtzCount(); i++)
	{
		Vertex *v = &_vertices[i];
		vec3 pos = make_vec3(v->pos[0], v->pos[1], v->pos[2]);
		vec3 t = mat4x4_mult_vec(model, pos);
		
		MtzFile *file = _ave->getMtz(i)->getMtzFile();
		if (add == 0)
		{
			file->setSelected(false);
		}
		
		if (file->isMarked())
		{
			continue;
		}

		if (t.x > x1 && t.x < x2 && t.y > y1 && t.y < y2)
		{
			file->setSelected(add >= 0 ? true : false);
		}
	}

	emit updateSelection();
	recolour();
}

void GLPoint::recolour()
{
	for (size_t i = 0; i < _ave->mtzCount(); i++)
	{
		Vertex &v = _vertices[i];
		v.color[0] = 0;
		v.color[1] = 0;
		v.color[2] = 0;

		MtzFile *file = _ave->getMtz(i)->getMtzFile();
		if (file->isSelected())
		{
			v.color[0] = 125;
			v.color[1] = 125;
		}
		else if (file->isMarked())
		{
			v.color[0] = 250;
			v.color[1] = 0;
		}
	}

	_keeper->update();
}

void GLPoint::repopulate()
{
	_vertices.clear();
	_indices.clear();

	for (size_t i = 0; i < _ave->mtzCount(); i++)
	{
		vec3 point = _ave->getPoint(i, _a, _b, _c);
		addPoint(point);
	}

	recolour();
	_keeper->update();
}

void GLPoint::setAverager(Averager *ave)
{
	_ave = ave;
	repopulate();
}

