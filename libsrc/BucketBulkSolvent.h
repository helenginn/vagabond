//
//  BucketBulkSolvent.h
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__BucketBulkSolvent__
#define __vagabond__BucketBulkSolvent__

#include <stdio.h>
#include "Bucket.h"

class BucketBulkSolvent : public Bucket
{
public:
	virtual void addSolvent();
	
	virtual ~BucketBulkSolvent() {}
private:
	void removeSlivers();
	bool sliverRemovalIteration(vec3 limits);
	void clearSliver(long i, long j, long k,
	                 long p, long q, long r);
};

#endif /* defined(__vagabond__BucketBulkSolvent__) */
