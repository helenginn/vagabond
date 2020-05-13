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

#include <string>
#include "GLPoint.h"
#include "Group.h"
#include <libsrc/FileReader.h>

GLPoint::GLPoint() : Plot3D()
{

}

void GLPoint::populate()
{
	for (size_t i = 0; i < _ave->mtzCount(); i++)
	{
		vec3 point = _ave->getPoint(i, _a, _b, _c);
		addPoint(point);
	}
}

size_t GLPoint::axisCount()
{
	return _ave->mtzCount();
}

std::string GLPoint::axisLabel(int i)
{
	double w = _ave->getDiagW(i);

	std::string str = f_to_str(w, 1);
	return str;
}
