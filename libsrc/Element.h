//
//  Element.h
//  vagabond
//
//  Created by Helen Ginn on 16/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Element__
#define __vagabond__Element__

#include <stdio.h>
#include <string>
#include <vector>
#include "shared_ptrs.h"

class Element
{
public:
	static void setupElements();

	Element(std::string symbol, std::string name, double covRadius, double mass);
	static ElementPtr getElement(std::string symbol);
	FFTPtr getDistribution();

	std::string getSymbol()
	{
		return _symbol;
	}

	std::string getName()
	{
		return _name;
	}

	double getCovalentRadius()
	{
		return _covRadius;
	}


private:
	std::string _symbol;
	std::string _name;
	double _covRadius;
	double _mass;
	double _density;
	FFTPtr _shape;

	static std::vector<ElementPtr> elements;

	void setCovalentRadius(double radius)
	{
		_covRadius = radius;
	}

	void setSymbol(std::string symbol)
	{
		_symbol = symbol;
	}

	void setName(std::string name)
	{
		_name = name;
	}

};

#endif /* defined(__vagabond__Element__) */
