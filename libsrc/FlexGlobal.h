//
//  FlexGlobal.hpp
//  vagabond
//
//  Created by Helen Ginn on 27/12/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#ifndef FlexGlobal_hpp
#define FlexGlobal_hpp

#include <stdio.h>
#include "shared_ptrs.h"

typedef enum
{
    FlexTargetMaximiseIsotropy,
    FlexTargetMatchOrigBFactor,
} FlexTarget;

class FlexGlobal
{
public:
    FlexGlobal();

    void setAtomGroup(AtomGroupPtr group)
    {
        _atomGroup = group;
    }

    static double score(void *object);

    void maximiseIsotropy();

    void matchOriginalBees()
    {
        _targetType = FlexTargetMatchOrigBFactor;
    }

    void setTargetBFactor(double value)
    {
        _targetIsoB = value;
    }
private:
    double _targetIsoB;

    double matchOriginalBeeScore();
    double maximiseIsotropyScore();
    AtomGroupPtr _atomGroup;

    FlexTarget _targetType;
};

#endif /* FlexGlobal_hpp */
