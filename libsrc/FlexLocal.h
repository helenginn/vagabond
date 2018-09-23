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
typedef std::map<BondPtr, double> BondCorrel;
typedef std::map<BondPtr, BondCorrel> BondBondCC;

typedef struct
{
	int index;
	int degree;
} BondDegree;

class FlexLocal
{
public:
	FlexLocal();
	
	void setPolymer(PolymerPtr pol, double shift);

	static double getScore(void *object);
	
	/** Performs preliminary work and refinement simultaneously */
	void refine();
	
	double getShift()
	{
		return _shift;
	}
private:
	void createAtomTargets();
	AtomTarget currentAtomValues();
	void createClustering();
	void reorganiseBondOrder();
	double bondRelationship(BondPtr bi, BondPtr bj);
	void scanBondParams();
	
	double directSimilarity();

	void chooseBestDifferenceThreshold();
	std::map<int, int> getClusterMembership(double threshold);

	PolymerPtr _polymer;

	AtomTarget _atomTargets;
	AtomTarget _atomOriginal;
	BondEffects _bondEffects;
	
	std::vector<AtomPtr> _atoms;
	std::vector<BondPtr> _bonds;
	std::vector<int> _reorderedBonds;
	std::vector<double> _b2bDiffs;
	std::vector<BondDegree> _degrees;
	std::map<BondPtr, int> _bondClusterIds;
	std::vector<ParamBandPtr> _paramBands;
	BondBondCC _bbCCs;
	
	int _afterBond;
	double _threshold;
	double _increment;
	double _anchorB;
	int _direct;
	int _window;
	int _run;
	double _shift;

};

#endif


