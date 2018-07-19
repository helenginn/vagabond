//
//  Chromosomal.h
//  vagabond
//
//  Created by Helen Ginn on 01/07/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#ifndef __vagabond__Chromosomal__
#define __vagabond__Chromosomal__

/* Abstract class for implementing genetic algorithms */

#include "shared_ptrs.h"

class Chromosomal
{
public:
	void evolve();
protected:
 	Chromosomal();

	virtual void randomise(double = 1.1) = 0;
	virtual void mutate() = 0;
	virtual double evaluate() = 0;
	virtual ChromosomalPtr makeCopy() = 0;
	
	virtual void haveSexWith(Chromosomal *other) = 0;
	virtual void copyOver(Chromosomal *_other) = 0;
	virtual void geneticCode() {};
	
private:
	int _sampled;

	void testPopulation();
	void breed();
	
	std::vector<ChromosomalPtr> _copies;
	std::vector<ChromosomalPtr> _survivors;
};

#endif
