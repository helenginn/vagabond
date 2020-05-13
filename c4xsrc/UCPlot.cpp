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

#include <string>
#include "UCPlot.h"
#include "Group.h"
#include "MtzFFT.h"
#include "MtzFile.h"
#include <libsrc/FileReader.h>

UCPlot::UCPlot() : Plot3D()
{

}

void UCPlot::populate()
{
	std::vector<double> ucs = _ave->getUnitCell();
	double mrw, mrf, srw, srf;
	_ave->averageRs(&mrw, &mrf, &srw, &srf);

	for (size_t i = 0; i < _ave->mtzCount(); i++)
	{
		std::vector<double> uc = _ave->getMtz(i)->getUnitCell();
		double rwork = _ave->getMtzFile(i)->getRWork() - mrw;
		double rfree = _ave->getMtzFile(i)->getRFree() - mrf;
		
		for (int j = 0; j < 6; j++)
		{
			uc[j] -= ucs[j];
		}
		
		uc.push_back(rwork);
		uc.push_back(rfree);

		vec3 point = make_vec3(uc[_a], uc[_b], uc[_c]);
		addPoint(point);
	}
}

size_t UCPlot::axisCount()
{
	return 8;
}

std::string UCPlot::axisLabel(int i)
{
	switch (i)
	{
		case 0:
		return "a";

		case 1:
		return "b";

		case 2:
		return "c";

		case 3:
		return "alpha";

		case 4:
		return "beta";

		case 5:
		return "gamma";
		
		case 6:
		return "rwork";

		case 7:
		return "rfree";
		
		default:
		return "?";
	}
}
