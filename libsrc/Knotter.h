//
//  Knotter.h
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Knotter__
#define __vagabond__Knotter__

#include <stdio.h>
#include "shared_ptrs.h"

/* Takes a side chain and ties all the atoms up together
 * with bonds! */

class Knotter
{
public:
	void setSidechain(SidechainPtr sidechain)
	{
		_sidechain = sidechain;
	}

	void setBackbone(BackbonePtr backbone)
	{
		_backbone = backbone;
	}

	void tie();
	void tieTowardsCTerminus();
private:
	SidechainPtr _sidechain;
	BackbonePtr _backbone;

	void makeGlycine();
	void makeCysteine();
	void makeValine();
	void makeSerine();
	void makeLysine();
	void makeThreonine();
	void makeHistidine();
	void makePhenylalanine();
	void makeMethionine();
};


#endif /* defined(__vagabond__Knotter__) */
