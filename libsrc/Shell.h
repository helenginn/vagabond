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

#ifndef __vagabond__Shell__
#define __vagabond__Shell__

#include <hcsrc/maths.h>

typedef struct
{
	double minRes;
	double maxRes;
	double rFactor;
	double scale;
	double std_err;
	double aveFo;
	double phi_spread;
	double count;
	std::vector<double> work1;
	std::vector<double> work2;
	std::vector<double> free1;
	std::vector<double> free2;
} ShellInfo;

inline ShellInfo makeShellInfo(double min, double max)
{
	ShellInfo shell;
	shell.minRes = min;
	shell.maxRes = max;
	shell.rFactor = 0;
	shell.scale = 1;
	shell.std_err = 0;
	shell.aveFo = 0;
	shell.phi_spread = 0;
	shell.count = 0;
	
	return shell;
}

/* resolutions are in real space */
inline int findShell(std::vector<ShellInfo> &shells, double res)
{
	if (res < shells.back().maxRes || 
	    res > shells.front().minRes)
	{
		return -1;
	}
	
	int size = shells.size();
	int min = 0;
	int max = size - 1;
	
	if (res >= shells[min].maxRes) return min;
	if (res <= shells[max].minRes) return max;
	
	while (true)
	{
		int chop = (max + min) / 2;
		int higher = (res <= shells[chop].maxRes);
		if (higher) 
		{
			min = chop;
		}
		else 
		{
			max = chop;
		}
		
		if (max - min == 1)
		{
			return max;
		}
	}
}

inline void makeShells(std::vector<ShellInfo> *shells, double min, double max,
                       int number = 20)
{
	/* Then apply to individual resolution bins */
	std::vector<double> bins;
	generateResolutionBins(min, max, number, &bins);

	/* Extend the final bin by a little bit, so as not to lose any
	 * stragglers. */
	bins[bins.size() - 1] *= 0.95;

	/* Make the series of shells */
	shells->clear();

	for (size_t i = 0; i < bins.size() - 1; i++)
	{
		ShellInfo shell = makeShellInfo(bins[i], bins[i + 1]);
		shells->push_back(shell);
	}

}

#endif
