//
//  Kabsch.hpp
//  vagabond
//
//  Created by Helen Ginn on 18/10/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#ifndef Kabsch_hpp
#define Kabsch_hpp

#include <stdio.h>
#include <vector>
#include "vec3.h"
#include "shared_ptrs.h"
#include "mat3x3.h"

class Kabsch
{
public:
	void setAtoms(std::vector<vec3> aAtoms, std::vector<vec3> bAtoms);
	mat3x3 run();
private:
	std::vector<std::vector<double> > _positions[2];
};

#endif /* Kabsch_hpp */
