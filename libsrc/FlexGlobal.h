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
#include "MapScoreWorkspace.h"
#include "RefinementStrategy.h"

/**
 * \class FlexGlobal
 * \brief Looks after scoring function for whole molecule movements.
 */

class FlexGlobal
{
public:
	FlexGlobal();

	/** Set atoms to consider in the target function. */
	void setAtomGroup(AtomGroupPtr group)
	{
		_atomGroup = group;
	}

	static double score(void *object);
	
	/** Crystal against which electron density should be matched. */
	void setCrystal(CrystalPtr crystal)
	{
		_crystal = crystal;	
	}
	
	void plot(std::string filename);
	
	double localCompareParams(Parameter &p1, Parameter &p2);
	double localScaleParam(Parameter &p1);
	static double scaleParam(void *obj, Parameter &p1);
	static double compareParams(void *obj, Parameter &p1, Parameter &p2);
private:
	bool _prepared;
	void prepareWorkspace();

	AtomGroupPtr _atomGroup;
	CrystalPtr _crystal;

	MapScoreWorkspace _workspace;
};

#endif /* FlexGlobal_hpp */
