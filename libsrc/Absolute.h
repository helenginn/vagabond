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

class Absolute : public Model, public std::enable_shared_from_this<Absolute>
{
public:
	Absolute(vec3 pos, double bFac);

	void addToMolecule(MoleculePtr molecule);
private:
	AtomPtr _atom;

	vec3 position;
	double bFactor;
};

#endif /* defined(__vagabond__Absolute__) */
