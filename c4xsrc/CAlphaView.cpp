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

#include "CAlphaView.h"
#include "MtzFile.h"
#include "Averager.h"
#include "KeeperGL.h"
#include <iostream>

CAlphaView::CAlphaView(MtzFile *mtz, vec3 centre)
{
	_mtzs.push_back(mtz);
	_centre = centre;
	_renderType = GL_LINES;
}

CAlphaView::CAlphaView(Averager *ave, vec3 centre)
{
	for (size_t i = 0; i < ave->mtzCount(); i++)
	{
		MtzFile *file = ave->getMtzFile(i);
		_mtzs.push_back(file);
	}
	_centre = centre;
	_renderType = GL_LINES;
}

void CAlphaView::initialisePrograms()
{
	GLObject::initialisePrograms(NULL, NULL);
}

void CAlphaView::repopulate()
{
	_indices.clear();
	_vertices.clear();
	_ends.clear();
	_starts.clear();
	vec3 nanVec = make_vec3(std::nan(""), std::nan(""), std::nan(""));

	for (size_t j = 0; j < _mtzs.size(); j++)
	{
		size_t begin = _vertices.size();
		std::vector<vec3> pos = _mtzs[j]->getAtomPositions();

		for (size_t i = 0; i < pos.size(); i++)
		{
			vec3 move = vec3_subtract_vec3(pos[i], _centre);
			addCAlpha(move);
		}
		
		addCAlpha(nanVec);

		size_t posSize = _vertices.size();
		
		_starts[_mtzs[j]] = begin;
		_ends[_mtzs[j]] = posSize;
	}
	
	_indices.pop_back();
	recolour();
}

void CAlphaView::recolour()
{
	for (size_t j = 0; j < _mtzs.size(); j++)
	{
		size_t begin = _starts[_mtzs[j]];
		size_t end = _ends[_mtzs[j]];
		
		for (size_t i = begin; i < end; i++)
		{
			_mtzs[j]->recolourVertex(&_vertices[i], true);
		}
	}
	
	_keeper->update();
}


void CAlphaView::addCAlpha(vec3 point)
{
	bool lastOK = (_vertices.size() > 0);
	
	if (lastOK)
	{
		Vertex last = _vertices.back();
		lastOK = (last.pos[0] == last.pos[0]);
	}
	
	if (point.x == point.x && lastOK)
	{
		_indices.push_back(_vertices.size() - 1);
		_indices.push_back(_vertices.size());
	}

	Vertex v;
	memset(v.pos, 0, sizeof(Vertex));
	
	v.color[3] = 1;
	v.pos[0] = point.x;
	v.pos[1] = point.y;
	v.pos[2] = point.z;

	_vertices.push_back(v);
}
