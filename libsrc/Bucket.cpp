//
//  Bucket.cpp
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Bucket.h"
#include <iostream>
#include "fftw3d.h"


void Bucket::scaleSolvent()
{
	if (!_data)
	{
		shout_at_helen("Need diffraction data to scale solvent");
	}
}

