//
//  maths.h
//  vagabond
//
//  Created by Helen Ginn on 22/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__vagmaths__
#define __vagabond__vagmaths__

#include <stdio.h>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <float.h>
#include "MapScoreWorkspace.h"

double correlation(std::vector<CoordVal> &vals);
double weighted_r_factor(std::vector<CoordVal> &vals);

#endif /* defined(__vagabond__maths__) */

