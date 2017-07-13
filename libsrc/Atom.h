//
//  Atom.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Atom__
#define __vagabond__Atom__

#include <stdio.h>
#include "shared_ptrs.h"
#include <vector>
#include "vec3.h"
#include "mat3x3.h"

class Atom
{
public:
	void addConnection(ModelPtr model);

	vec3 getPosition()
	{
		return _position;
	}

	void setPosition(vec3 pos)
	{
		_position = pos;
	}

	void addToMap(FFTPtr fft, mat3x3 unit_cell);

private:
	std::vector<ModelPtr> connections;

	vec3 _position;
};

#endif /* defined(__vagabond__Atom__) */
