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

#include "AveUnitCell.h"
#include "MtzFFT.h"
#include "MtzFile.h"
#include <libsrc/maths.h>

AveUnitCell::AveUnitCell(Group *group) : Average(group)
{

}

void AveUnitCell::calculate()
{
	_ucs = std::vector<double>(6, 0.);
	int count = 0;

	for (size_t i = 0; i < _mtzs.size(); i++)
	{
		if (_mtzs[i]->getMtzFile()->isDead())
		{
			continue;
		}

		std::vector<double> uc = _mtzs[i]->getUnitCell();
		
		if (uc.size() < 6)
		{
			continue;
		}
		
		for (int j = 0; j < 6; j++)
		{
			_ucs[j] += uc[j];
		}
		
		count++;
	}
	
	if (count == 0)
	{
		count = 1;
	}

	for (int j = 0; j < 6; j++)
	{
		_ucs[j] /= (double)count;
	}
}

double AveUnitCell::findCorrelation(MtzFFTPtr one, MtzFFTPtr two)
{
	std::vector<double> uc1 = one->getUnitCell();
	std::vector<double> uc2 = two->getUnitCell();
	
	CorrelData cd = empty_CD();

	for (int j = 0; j < 6; j++)
	{
		uc1[j] -= _ucs[j];
		uc2[j] -= _ucs[j];
		
		add_to_CD(&cd, uc1[j], uc2[j]);
	}
	
	double correl = evaluate_CD(cd);
	return correl;
}
