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

#include "Average.h"
#include "Group.h"
#include <iostream>

Average::Average(Group *group)
{
	_group = group;
	_symmetric = true;
	
	if (_group)
	{
		_mtzs = group->mtzs();
	}
}

void Average::findIntercorrelations(Group *other, double **svd)
{
	if (_symmetric)
	{
		for (size_t i = 1; i < other->mtzCount(); i++)
		{
			for (size_t j = 0; j < i; j++)
			{
				double cc = 0;

				MtzFFTPtr one = other->getMtz(i);
				MtzFFTPtr two = other->getMtz(j);

				cc = findCorrelation(one, two);

				svd[i][j] = cc;
				svd[j][i] = cc;
			}

			std::cout << "." << std::flush;
		}
	}
	else
	{
		for (size_t i = 0; i < other->mtzCount(); i++)
		{
			for (size_t j = 0; j < other->mtzCount(); j++)
			{
				if (i == j)
				{
					continue;
				}

				MtzFFTPtr one = other->getMtz(i);
				MtzFFTPtr two = other->getMtz(j);

				double cc1 = findCorrelation(one, two);
				double cc2 = findCorrelation(two, one);
				
				if (cc1 != cc1 && cc2 != cc2)
				{
					svd[i][j] = nan(" ");
				}
				else if (fabs(cc1) < 1e-6 || cc1 != cc1)
				{
					svd[i][j] = cc2;
				}
				else if (fabs(cc2) < 1e-6 || cc2 != cc2)
				{
					svd[i][j] = cc1;
				}
				else
				{
					svd[i][j] = (cc1 + cc2) / 2;
				}
			}

			std::cout << "." << std::flush;
		}
	}

	std::cout << std::endl;

}
