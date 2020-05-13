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

#include "PlotR.h"
#include "MtzFile.h"
#include "Group.h"

PlotR::PlotR() : Plot3D()
{

}

void PlotR::populate()
{
	for (size_t i = 0; i < _ave->mtzCount(); i++)
	{
		MtzFile *f = _ave->getMtzFile(i);
		vec3 point = empty_vec3();
		point.x = f->getRWork();
		point.y = f->getRFree();

		addPoint(point);
	}
}

size_t PlotR::axisCount()
{
	return 2;
}

std::string PlotR::axisLabel(int i)
{
	switch (i)
	{
		case 0:
		return "r_work";

		case 1:
		return "r_free";

		case 2:
		return "[none]";
		
		default:
		return "?";
	}
}
