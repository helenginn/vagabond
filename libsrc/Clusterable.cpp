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
#include "RefinementStrategy.h"

Clusterable::Clusterable(const int ndims)
{
	_ndims = ndims;
	
	for (int n = 0; n < _ndims; n++)
	{
		RefinableDoublePtr object;
		object = RefinableDoublePtr(new RefinableDouble(n, this));
		double value = random() / (double)RAND_MAX;
		RefinableDouble::setDouble(&*object, value);
		object->setGradientFunction(gradientFunc);
		
		_coords.push_back(object);
	}
}

double Clusterable::dotWith(ClusterablePtr cluster)
{
	double sum = 0;

	for (int i = 0; i < _ndims; i++)
	{
		double contrib = RefinableDouble::getDouble(&cluster->_coords[i]);
		contrib *= RefinableDouble::getDouble(&_coords[i]);
		
		sum += contrib;
	}
	
	return sum;
}

double Clusterable::sumContributionToEval()
{
	double fx = 0;
	
	for (int i = 0; i < _leftClusters.size(); i++)
	{
		ClusterablePtr cluster = _leftClusters[i].lock();
		double dot = dotWith(cluster);
		double cc = _ccMap[cluster];
		
		fx += (cc - dot);
	}
	
	return fx;
}

double Clusterable::gradient(int tag)
{

}

double Clusterable::gradientFunc(void *object, int tag)
{
	Clusterable *clust = static_cast<Clusterable *>(object);
	double grad = clust->gradient(tag);
	return grad;
}

void Clusterable::addRelationship(ClusterablePtr cluster, double cc)
{
	ClusterablePtr me = shared_from_this();
	_ccMap[cluster] = cc;
	cluster->_ccMap[me] = cc;
	cluster->_rightClusters.push_back(me);
	
	_leftClusters.push_back(cluster);
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

