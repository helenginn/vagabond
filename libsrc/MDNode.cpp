//
//  MDNode.cpp
//  vagabond
//
//  Created by Helen Ginn on 04/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "Shouter.h"
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
	if (_nodes == NULL)
	{
		_value += value;
		_weight += 1;
		return;
	}
	
	int node = 0;

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
			return;		
		}

		if (side == 1)
		{
			node += pow(2, i);
		}
	}
	
	if (node >= 0)
	{
		_nodes[node]->addToNode(dimvals, value);
	}
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
	int total = pow(2, _dims);
	
	for (int i = 0; i < total; i++)
	{
		delete _nodes[i];
	}
	
	delete [] _nodes;
	
	delete [] _mins;
	delete [] _maxes;
}
