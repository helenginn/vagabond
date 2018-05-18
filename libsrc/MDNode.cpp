//
//  MDNode.cpp
//  vagabond
//
//  Created by Helen Ginn on 04/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "Shouter.h"
#include "Plucker.h"
#include "MDNode.h"
#include <iostream>
#include <vector>
#include "FileReader.h"
#include "CSV.h"
#include <math.h>

MDNode::MDNode(int dims)
{
	if (dims < 0)
	{
		shout_at_helen("Wtf");
	}
	
	_dims = dims;
	_mins = new double[dims];
	_maxes = new double[dims];
	
	_value = 0;
	_weight = 0;
	
	_plucker = NULL;
	_nodes = NULL;
}

void MDNode::setDimension(int dim, double min, double max)
{
	if (dim >= _dims)
	{
		return;
	}
	
	_mins[dim] = min;
	_maxes[dim] = max;
}

void MDNode::splitNode(int divisions, int remaining)
{
	if (remaining <= 0)
	{
		return;
	}

	/* This many new nodes are required */
	int total = pow(2, _dims); /* 1 << divisions?? */
	_nodes = new MDNode *[total];
	
	/* Each node will occupy the left or right half of each dimension.
	* We need to establish which sides we should be taking for each one. */
	for (int i = 0; i < total; i++)
	{
		_nodes[i] = new MDNode(_dims);

		for (int j = 0; j < _dims; j++)
		{
			int test = pow(2, j);
			
			int side = (test & i) > 0;
			
			double min = _mins[j];
			double max = _maxes[j];
			double ave = (_maxes[j] + _mins[j]) / 2;
			
			if (side == 1)
			{
				min = ave;
			}
			else
			{
				max = ave;
			}
			
			_nodes[i]->setDimension(j, min, max);
		}
		
		_nodes[i]->splitNode(divisions, remaining - 1);
	}
}

void MDNode::addToNode(double *dimvals, double value)
{
	MDNode *node = findNode(dimvals);
	
	if (node)
	{
		node->_value += value;
		node->_weight += 1;
	}
}

MDNode *MDNode::findNode(double *dimvals)
{
	int node = 0;
	
	if (_nodes == NULL)
	{
		return this;
	}

	for (int i = 0; i < _dims; i++)
	{
		double min = _mins[i];
		double max = _maxes[i];
		
		double ave = (_maxes[i] + _mins[i]) / 2;
		int side = -1;
		
		if (dimvals[i] > _mins[i] && dimvals[i] < ave)
		{
			side = 0;		
		}
		else if (dimvals[i] >= ave && dimvals[i] < _maxes[i])
		{
			side = 1;	
		}
		
		if (side == -1)
		{
			return NULL;		
		}

		if (side == 1)
		{
			node += pow(2, i);
		}
	}
	
	return _nodes[node]->findNode(dimvals);
}

void MDNode::addNodesToPlucker(Plucker *pluck, int dim,
                               double value, double subtract)
{
	double min = _mins[dim];
	double max = _maxes[dim];

	if (value < min || value > max)
	{
		return;
	}

	for (int i = 0; i < nodeCount(); i++)
	{
		node(i)->addNodesToPlucker(pluck, dim, value, subtract);
	}
	
	if (_nodes != NULL)
	{
		return;
	}

	double val = getValue();
	val -= subtract;

	pluck->addPluckable(this, val);
}

void MDNode::makePlucker(int dim, double value, double subtract)
{
	if (_plucker)
	{
		delete _plucker;
		_plucker = NULL;
	}
	
	_plucker = new Plucker();
	_plucker->setGranularity(0.04);

	addNodesToPlucker(_plucker, dim, value, 0.20);
}

void MDNode::addChildrenToCSV(CSV *csv)
{
	if (_nodes != NULL)
	{
		int total = pow(2, _dims);
	
		for (int i = 0; i < total; i++)
		{
			_nodes[i]->addChildrenToCSV(csv);
		}
	}
	else
	{
		std::vector<double> entry;
		
		for (int i = 0; i < _dims; i++)
		{
			double ave = (_mins[i] + _maxes[i]) / 2;
			entry.push_back(ave);
		}
		
		entry.push_back(_value / _weight);
		
		csv->addEntry(entry);
	}
}

void MDNode::addToCSV(CSV *csv)
{
	for (int i = 0; i < _dims; i++)
	{
		csv->addHeader("dim" + i_to_str(i));
	}
	
	csv->addHeader("value");
	
	addChildrenToCSV(csv);
}

MDNode::~MDNode()
{
	if (_plucker != NULL)
	{
		delete _plucker;
		_plucker = NULL;
	}
	
	if (_nodes != NULL)
	{
		int total = pow(2, _dims);

		for (int i = 0; i < total; i++)
		{
			delete _nodes[i];
		}

		delete [] _nodes;
	}
	
	delete [] _mins;
	delete [] _maxes;
}
