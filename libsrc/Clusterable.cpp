//
//    Clusterable.cpp
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

#include "Clusterable.h"
#include "FileReader.h"
#include "RefinableDouble.h"
#include "RefinementStrategy.h"

Clusterable::Clusterable(const int ndims)
{
	_ndims = ndims;
	
	for (int n = 0; n < _ndims; n++)
	{
		RefinableDouble object = RefinableDouble(n, this);
		double value = random() / (double)RAND_MAX;
		RefinableDouble::setDouble(&object, value);
		
		_coords.push_back(object);
	}
}

void Clusterable::addRelationship(ClusterablePtr cluster, double cc)
{
	ClusterablePtr me = shared_from_this();
	_ccMap[cluster] = cc;
	cluster->_ccMap[me] = cc;
}

void Clusterable::addParamsToStrategy(RefinementStrategyPtr strategy)
{
	for (int n = 0; n < _ndims; n++)
	{
		strategy->addParameter(&_coords[n], RefinableDouble::getDouble,
		                       RefinableDouble::setDouble,
		                       0.5, 0.01, "coord_" + i_to_str(n),
		                       RefinableDouble::getGradient);
	}
}

