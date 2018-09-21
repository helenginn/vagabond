//
//    Clusterable.h
//    Vagabond - stuff

//    Copyright (C) 2018  Helen Ginn
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __vagabond__Clusterable__
#define __vagabond__Clusterable__

#include "shared_ptrs.h"
#include "RefinableDouble.h"
#include <map>

class RefinableDouble;

class Clusterable : public boost::enable_shared_from_this<Clusterable>
{
public:
	Clusterable(const int ndims);
	void addParamsToStrategy(RefinementStrategyPtr strategy);
	
	/* Only add in the case that i < j, the appropriate other
	 * relationship will be set up automatically. */
	void addRelationship(ClusterablePtr cluster, double cc);
	
	static double gradientFunc(void *object, int tag);
	double sumContributionToEval();
	
	double getCoord(int i)
	{
		return RefinableDouble::getDouble(&*_coords[i]);
	}
	
	double ccWith(ClusterablePtr cluster);
	double dotWith(ClusterablePtr cluster);
	
	std::vector<double> coords();
private:
	int _ndims;
	std::vector<RefinableDoublePtr> _coords;

	double gradient(int tag);

	std::map<ClusterableWkr, double> _ccMap;
	std::vector<ClusterableWkr> _leftClusters, _rightClusters;
};

#endif

