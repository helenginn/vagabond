//
//  Absolute.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Absolute__
#define __vagabond__Absolute__

#include <stdio.h>
#include "Model.h"
#include "vec3.h"
#include <string>

class Absolute : public Model
{
public:
	Absolute(vec3 pos, double bFac, std::string element, double occValue);


// Model virtual functions:
	virtual FFTPtr getDistribution();
	virtual void addToMolecule(MoleculePtr molecule);

private:
	AtomPtr _atom;
	std::string _element;
	double _occupancy;

	vec3 position;
	double bFactor;
};

#endif /* defined(__vagabond__Absolute__) */
