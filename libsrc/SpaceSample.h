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
#include <hcsrc/vec3.h>
#include "shared_ptrs.h"
#include <vector>

class ConfSpace;
class Superpose;

typedef std::vector<double> SpacePoint;

class SpaceSample
{
public:
	SpaceSample(ConfSpace *space);
	
	void addDiagonalParameters(RefinementStrategyPtr strategy);
	void addTensorParameters(RefinementStrategyPtr strategy);
	void addAverageParameters(RefinementStrategyPtr strategy);
	void generatePoints(TotalModelPtr total, bool use = true);
	double getDeviation(int resi, int sample, bool phi);

	void savePositions();
	
	void calculateSuperpositions();
	
	bool hasPoints()
	{
		return _points.size() > 0;
	}
	
	vec3 point3D(int i);
	
	~SpaceSample();
	
	void setAtoms(AtomList atoms);
	
	Superpose *superpose()
	{
		return _superpose;
	}

	size_t motionCount()
	{
		return _motions.size();
	}
	
	MotionPtr getMotion(int i)
	{
		return _motions[i];
	}
	
	void addMotion(AnchorPtr a, MotionPtr mot)
	{
		_anchors.push_back(a);
		_motions.push_back(mot);
	}
private:
	void fillInTensorGaps();
	ConfSpace *_space;
	int _used;

	double *_average;
	double *_stdev;
	HelenCore::Matrix _tensor;
	Superpose *_superpose;
	
	std::vector<MotionPtr> _motions;
	std::vector<AnchorPtr> _anchors;
	std::vector<SpacePoint> _points;
	std::vector<SpacePoint> _original;
	std::vector<AnyPtr> _anys;
};

#endif
