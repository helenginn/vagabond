// vagabond
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

#ifndef __vagabond__SpaceSample__
#define __vagabond__SpaceSample__

#include <hcsrc/Matrix.h>
#include "shared_ptrs.h"
#include <vector>

class ConfSpace;

typedef std::vector<double> SpacePoint;

class SpaceSample
{
public:
	SpaceSample(ConfSpace *space);
	
	void addDiagonalParameters(RefinementStrategyPtr strategy);
	void addTensorParameters(RefinementStrategyPtr strategy);
	void addAverageParameters(RefinementStrategyPtr strategy);
	void generatePoints(CrystalPtr cryst, int m = 0);
	double getTorsionDeviation(int resi, int sample);
	double getWhackDeviation(int resi, int sample);
	
	bool hasPoints()
	{
		return _points.size() > 0;
	}
	
	~SpaceSample();

private:
	ConfSpace *_space;

	double *_average;
	double *_stdev;
	HelenCore::Matrix _tensor;
	
	std::vector<SpacePoint> _points;
	std::vector<AnyPtr> _anys;
};

#endif
