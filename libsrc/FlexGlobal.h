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
	FlexTargetMatchElectronDensity,
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

	void matchElectronDensity()
	{
		_targetType = FlexTargetMatchElectronDensity;
	}

	void matchOriginalBees()
	{
		_targetType = FlexTargetMatchOrigBFactor;
	}

	void setTargetBFactor(double value)
	{
		_targetIsoB = value;
	}
	
	void setCrystal(CrystalPtr crystal)
	{
		_crystal = crystal;	
	}
private:
	double _targetIsoB;

	double matchOriginalBeeScore();
	double maximiseIsotropyScore();
	double matchElectronDensityScore();
	AtomGroupPtr _atomGroup;
	CrystalPtr _crystal;

	FlexTarget _targetType;
};

#endif /* FlexGlobal_hpp */
