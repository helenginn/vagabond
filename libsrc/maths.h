//
//  maths.h
//  vagabond
//
//  Created by Helen Ginn on 22/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__maths__
#define __vagabond__maths__

#include <stdio.h>
#include <vector>

double scale_factor(std::vector<double> &set1, std::vector<double> &set2);
double r_factor(std::vector<double> &set1, std::vector<double> &set2);
double correlation(std::vector<double> &vec1, std::vector<double> &vec2);

typedef double (*two_dataset_op)(std::vector<double>&, std::vector<double>&);

#endif /* defined(__vagabond__maths__) */
