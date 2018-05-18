//
//  BucketMonteCarlo.h
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__BucketMonteCarlo__
#define __vagabond__BucketMonteCarlo__

#include <stdio.h>
#include "Bucket.h"

class BucketMonteCarlo : public Bucket
{
public:
	BucketMonteCarlo();
	virtual void addSolvent();
	
	virtual ~BucketMonteCarlo() {}
private:
	void loadAnalysis(std::string filename);
	void makePluckers(double distance);

	AtomGroupPtr allWaters();
	AtomGroupPtr makeNewWaters(Plucker *plucker, int total);
	
	AtomGroupPtr _baseGrp;
	AtomGroupPtr _extra;
	
	int _loadedAnalysis;
};



#endif /* defined(__vagabond__BucketMonteCarlo) */
