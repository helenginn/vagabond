// Vagabond
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

#ifndef __vagabond__SimXW__
#define __vagabond__SimXW__

// relationship between X and W from Sim 1959

namespace Vagabond
{
	static int simXs[] = {0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0,
	10.0};
	static int simWs[] = {0, 0.157, 0.308, 0.443, 0.561, 0.738, 0.850,
	0.952, 0.982, 0.992, 1.0};
	static size_t xwSize = 11;
	
	double getW(double x)
	{
		double w = 1.;
		for (size_t i = 0; i < xwSize - 1; i++)
		{
			if (x >= simXs[i] && x < simXs[i + 1])
			{
				double width = (simXs[i + 1] - simXs[i]);
				double prop  = (x - simXs[i]) / width;
				
				w = simWs[i] + prop * (simWs[i + 1] - simWs[i]);
				break;
			}
		}

		return w;
	}
};

#endif
