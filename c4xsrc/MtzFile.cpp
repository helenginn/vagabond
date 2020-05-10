// clusterxxxx
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

#include "MtzFile.h"
#include <QFont>
#include "QuickAtoms.h"

MtzFile::MtzFile(std::string filename)
{
	_filename = filename;
	_mark = false;
	_sele = false;
	_dead = false;
	_quickAtoms = new QuickAtoms(this);
}

void MtzFile::recolourVertex(Vertex *v, bool fullDead)
{
	v->color[0] = 0;
	v->color[1] = 0;
	v->color[2] = 0;

	if (isDead())
	{
		v->color[0] = 100. / 255.;
		v->color[1] = 100. / 255.;
		v->color[2] = 100. / 255.;
		
		if (fullDead)
		{
			v->color[3] = 0.;
		}
	}
	if (isSelected())
	{
		v->color[0] = 200. / 255.;
		v->color[1] = 200. / 255.;
	}
	if (isMarked())
	{
		v->color[0] = 255 / 255;
		v->color[1] = 0;
	}
}
