//
//  NelderMead.h
//   cppxfel - a collection of processing algorithms for XFEL diffraction data.

//    Copyright (C) 2017  Helen Ginn
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __cppxfel__NelderMead__
#define __cppxfel__NelderMead__

#include <stdio.h>
#include "shared_ptrs.h"
#include "RefinementStrategy.h"

typedef std::pair<std::vector<double>, double> TestPoint;

class NelderMead : public RefinementStrategy
{
private:
    double alpha;
    double gamma;
    double rho;
    double sigma;
    
    std::vector<TestPoint> testPoints;
    
    void setWorstTestPoint(TestPoint &newPoint);
    TestPoint *worstTestPoint();
    void orderTestPoints();
    void evaluateTestPoint(int num);
    void evaluateTestPoint(TestPoint *testPoint);
    void setTestPointParameters(TestPoint *testPoint);
    std::vector<double> calculateCentroid();
    
    TestPoint reflectOrExpand(std::vector<double> centroid, double scale);
    TestPoint reflectedPoint(std::vector<double> centroid);
    TestPoint expandedPoint(std::vector<double> centroid);
    TestPoint contractedPoint(std::vector<double> centroid);
    void reduction();
    
    void addPoints(std::vector<double> *point, std::vector<double> pointToAdd);
    void scalePoint(std::vector<double> *point, double scale);
    void subtractPoints(std::vector<double> *point, std::vector<double> pointToSubtract);
public:
    void init();
    NelderMead() : RefinementStrategy() { init(); };
    virtual void refine();
    
    virtual void clearParameters();
};

#endif /* defined(__cppxfel__NelderMead__) */
