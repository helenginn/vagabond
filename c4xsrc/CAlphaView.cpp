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

struct _Vertex;
typedef _Vertex Vertex;

#include "CAlphaView.h"
#include "MtzFile.h"
#include "Group.h"
#include "KeeperGL.h"
#include "QuickAtoms.h"
#include <iostream>
#include <libsrc/FileReader.h>

CAlphaView::CAlphaView(MtzFile *mtz, vec3 centre)
{
	_mtzs.push_back(mtz);
	_centre = centre;
	_renderType = GL_LINES;
}

CAlphaView::CAlphaView(Group *ave)
{
	for (size_t i = 0; i < ave->mtzCount(); i++)
	{
		MtzFile *file = ave->getMtzFile(i);
		_mtzs.push_back(file);
	}

	_centre = ave->getCentre();
	_renderType = GL_LINES;
}

void CAlphaView::updateRs()
{
	CorrelData cd = empty_CD();

	for (size_t i = 0; i < _mtzs.size(); i++)
	{
		if (_mtzs[i]->isDead())
		{
			continue;
		}

		add_to_CD(&cd, _mtzs[i]->getRWork(), _mtzs[i]->getRFree());
	}

	means_stdevs_CD(cd, &_mean_rwork, &_mean_rfree, 
	                &_stdev_rwork, &_stdev_rfree);
}

void CAlphaView::initialisePrograms()
{
	SlipObject::initialisePrograms(NULL, NULL);
}

void CAlphaView::repopulate()
{
	_indices.clear();
	_vertices.clear();

	_ends.clear();
	_starts.clear();

	for (size_t j = 0; j < _mtzs.size(); j++)
	{
		size_t begin = _vertices.size();
		_mtzs[j]->getQuickAtoms()->populateCAlphaView(this);
		size_t posSize = _vertices.size();
		
		_starts[_mtzs[j]] = begin;
		_ends[_mtzs[j]] = posSize;
	}
	
	_indices.pop_back();
	recolour();
	updateRs();
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
	
	_c4xKeeper->update();
}


void CAlphaView::addCAlpha(vec3 point)
{
	bool lastOK = (_vertices.size() > 0);
	
	if (lastOK)
	{
		Helen3D::Vertex last = _vertices.back();
		lastOK = (last.pos[0] == last.pos[0]);
	}
	
	if (point.x == point.x && lastOK)
	{
		_indices.push_back(_vertices.size() - 1);
		_indices.push_back(_vertices.size());
	}

	Helen3D::Vertex v;
	memset(v.pos, 0, sizeof(GLfloat) * 3);
	vec3 move = vec3_subtract_vec3(point, _centre);
	
	v.color[3] = 1;
	v.pos[0] = move.x;
	v.pos[1] = move.y;
	v.pos[2] = move.z;

	_vertices.push_back(v);
}

std::string CAlphaView::getRworkRfree()
{
	std::string str;
	str += "Rwork: " + f_to_str(_mean_rwork * 100, 1);
	
	if (_mtzs.size() > 1)
	{
		str += " (stdev ";
		str += f_to_str(_stdev_rwork * 100, 1) + ")";
	}

	str += " %\nRfree: " + f_to_str(_mean_rfree * 100, 1);
	
	if (_mtzs.size() > 1)
	{
		str += " (stdev ";
		str += f_to_str(_stdev_rfree * 100, 1) + ")";
	}

	str += " %\n";
	return str;
}
