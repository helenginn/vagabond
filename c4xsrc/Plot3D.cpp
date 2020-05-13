// cluster4x
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

#include "KeeperGL.h"
#include "Plot3D.h"
#include "shaders/Blob_fsh.h"
#include "shaders/Blob_vsh.h"
#include "MtzFFT.h"
#include "MtzFile.h"
#include "Group.h"

Plot3D::Plot3D() : GLObject()
{
	_ave = NULL;
	_keeper = NULL;
	_renderType = GL_POINTS;
	_a = 0;
	_b = 1;
	_c = 2;
}

void Plot3D::initialisePrograms()
{
	std::string fsh = blobFsh();
	std::string vsh = blobVsh();
	GLObject::initialisePrograms(&vsh, &fsh);
}

void Plot3D::addPoint(vec3 point)
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

void Plot3D::selectInWindow(float x1, float y1, float x2, float y2,
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

void Plot3D::repopulate()
{
	_vertices.clear();
	_indices.clear();

	populate();

	recolour();
	_keeper->update();
}

void Plot3D::setGroup(Group *ave)
{
	_ave = ave;
	repopulate();
}


void Plot3D::recolour()
{
	for (size_t i = 0; i < _ave->mtzCount(); i++)
	{
		Vertex *v = &_vertices[i];
		MtzFile *file = _ave->getMtz(i)->getMtzFile();
		file->recolourVertex(v);
	}

	_keeper->update();
}

