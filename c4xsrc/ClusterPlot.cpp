// Cluster4x
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

#include "ClusterPlot.h"
#include "KeeperGL.h"
#include "Group.h"
#include "MtzFile.h"
#include "MtzFFT.h"
#include "shaders/Blob_fsh.h"
#include "shaders/Blob_vsh.h"

ClusterPlot::ClusterPlot() : Plot3D()
{
	_renderType = GL_POINTS;
	_grp = NULL;
	_keeper = NULL;
	setName("ClusterPlot");
	_fString = pointFsh();
	_vString = pointVsh();
}

void ClusterPlot::selectInWindow(float x1, float y1, float x2, float y2,
                             int add)
{
	for (size_t i = 0; i < _grp->mtzCount(); i++)
	{
		Helen3D::Vertex *v = &_vertices[i];
		vec3 pos = make_vec3(v->pos[0], v->pos[1], v->pos[2]);
		double last = 1;
		vec3 model = mat4x4_mult_vec3(_model, pos, &last);
		vec3 t = mat4x4_mult_vec3(_proj, model, &last);
		vec3_mult(&t, 1 / last);
		
		MtzFile *file = _grp->getMtz(i)->getMtzFile();

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
	
	Plot3D::selectInWindow(x1, y1, x2, y2, add);
}

void ClusterPlot::recolour()
{
	for (size_t i = 0; i < _grp->mtzCount(); i++)
	{
		Helen3D::Vertex *v = &_vertices[i];
		MtzFile *file = _grp->getMtz(i)->getMtzFile();
		file->recolourVertex(v);
		
		if (_texts.size() > i)
		{
			_texts[i]->setColour(v->color[0], v->color[1], 
			                     v->color[2], v->color[3]);
			_texts[i]->prepare();
		}
	}

	_keeper->update();
}

void ClusterPlot::extraUniforms()
{
	Plot3D::extraUniforms();
}
