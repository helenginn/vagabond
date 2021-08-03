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

bool Average::_vectorList = false;

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
	int total = 0;
	int missing = 0;

	std::cout << "find intercorrel Vectorlist: " << _vectorList << std::endl;
	
	if (_vectorList)
	{
		for (size_t i = 0; i < other->mtzCount(); i++)
		{
			for (size_t j = 0; j < other->mtzCount(); j++)
			{
				if (i == j)
				{
//					continue;
				}

				MtzFFTPtr one = other->getMtz(i);
				MtzFFTPtr two = other->getMtz(j);

				double cc = findCorrelation(one, two);

				svd[i][j] = cc;
			}
		}

		return;
	}

	if (_symmetric)
	{
		for (size_t i = 1; i < other->mtzCount(); i++)
		{
			for (size_t j = 0; j < i; j++)
			{
				double cc = 0;
				total += 2;

				MtzFFTPtr one = other->getMtz(i);
				MtzFFTPtr two = other->getMtz(j);

				cc = findCorrelation(one, two);
				
				if (cc != cc)
				{
					missing += 2;
				}

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
				total++;

				MtzFFTPtr one = other->getMtz(i);
				MtzFFTPtr two = other->getMtz(j);

				double cc1 = findCorrelation(one, two);
				double cc2 = findCorrelation(two, one);
				
				if (cc1 != cc1 && cc2 != cc2)
				{
					svd[i][j] = nan(" ");
					missing++;
				}
				else if (cc1 == cc1)
				{
					svd[i][j] = cc1;
				}
				else if (cc1 != cc1)
				{
					svd[i][j] = cc2;
				}
				else if (cc2 != cc2)
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

	if (missing > 0)
	{
		std::cout << "Total " << total << " matrix pairs, " <<
		missing << " missing entries." << std::endl;
	}

	std::cout << std::endl;

}
