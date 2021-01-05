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
#include <hcsrc/RefinementStrategy.h>

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
	
	void plot(std::string filename = "");
	
	void setReciprocalRefinement(bool recip = true)
	{
		_recip = recip;
	}
	
	MapScoreWorkspace &getWorkspace()
	{
		return _workspace;
	}
	
	void reportTimings();
	void recalculateConstant();
	void prepareComparisons(RefinementStrategyPtr str);
	static double compareParams(void *object, Parameter &p1, Parameter &p2);
private:
	bool _prepared;
	bool _recip;
	void prepareWorkspace();
	static int _plotCycle;

	AtomGroupPtr _atomGroup;
	CrystalPtr _crystal;

	MapScoreWorkspace _workspace;
	std::map<int, std::map<int, double> > _paramCorrel;
	std::vector<Parameter> _paramList;
	std::vector<VagFFTPtr> _mapList;
};

#endif /* FlexGlobal_hpp */
