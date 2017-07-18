//
//  BucketUniform.h
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__BucketUniform__
#define __vagabond__BucketUniform__

#include <stdio.h>
#include "Bucket.h"

class BucketUniform : public Bucket
{
public:
	virtual void addSolvent(FFTPtr map);
private:
};

#endif /* defined(__vagabond__BucketUniform__) */
