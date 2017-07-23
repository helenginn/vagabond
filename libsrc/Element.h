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
#include "ScatterFactors.h"

class Element
{
public:
	static void setupElements();

	Element(std::string symbol, std::string name);
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


private:
	std::string _symbol;
	std::string _name;
	float _scattering[ScatterFactors::numScatter];
	FFTPtr _shape;

	static std::vector<ElementPtr> elements;

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
