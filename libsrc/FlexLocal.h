//
//  FlexLocal.h
//  vagabond
//
//  Created by Helen Ginn, 2018
//  Copyright © 2018 Helen Ginn. All rights reserved.
//

#ifndef FlexLocal_h
#define FlexLocal_h

#include "shared_ptrs.h"
#include <map>

typedef std::map<AtomPtr, double> AtomTarget;
typedef std::map<BondPtr, AtomTarget> BondEffects;
typedef std::map<BondPtr, ClusterablePtr> BondClusters;

class FlexLocal
{
public:
	FlexLocal();
	
	void setPolymer(PolymerPtr pol)
	{
		_polymer = pol;
	}

	static double getScore(void *object);
	
	/** Performs preliminary work and refinement simultaneously */
	void refine();
	static double clusterScore(void *object);
private:
	void createAtomTargets();
	AtomTarget currentAtomValues();
	void createClustering();
	double bondRelationship(BondPtr bi, BondPtr bj);
	void scanBondParams();
	double clusterScore();
	
	double directSimilarity();

	PolymerPtr _polymer;

	AtomTarget _atomTargets;
	AtomTarget _atomOriginal;
	BondEffects _bondEffects;
	
	std::vector<AtomPtr> _atoms;
	std::vector<BondPtr> _bonds;
	BondClusters _clusters;
	
	double _anchorB;
	int _nDims;
	int _direct;
	int _window;
	int _run;
	double _shift;

};

#endif


