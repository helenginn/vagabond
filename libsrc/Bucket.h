//
//  Bucket.h
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Bucket__
#define __vagabond__Bucket__

#include <stdio.h>
#include "shared_ptrs.h"
#include <vector>
#include "vec3.h"

class Bucket
{
public:
	virtual void addSolvent(FFTPtr map) = 0;

protected:
	void findBulkSolvent(FFTPtr map);

private:
	void makeCheckingShifts(FFTPtr map);

	std::vector<vec3> fullCheckVecs;
};

#endif /* defined(__vagabond__Bucket__) */
