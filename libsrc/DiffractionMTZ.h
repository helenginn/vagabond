//
//  DiffractionMTZ.h
//  vagabond
//
//  Created by Helen Ginn on 22/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__DiffractionMTZ__
#define __vagabond__DiffractionMTZ__

#include <stdio.h>
#include "Diffraction.h"

class DiffractionMtz : public Diffraction
{
public:
	virtual void load();

	void syminfoCheck();
private:
	
};

#endif /* defined(__vagabond__DiffractionMTZ__) */
