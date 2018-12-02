//
//  FlexLocal.h
//  vagabond
//
//  Created by Helen Ginn, 2018
//  Copyright © 2018 Helen Ginn. All rights reserved.
//

#ifndef FlexLocal_h
#define FlexLocal_h

#include "RefinementStrategy.h"
#include "shared_ptrs.h"
#include <map>

typedef std::map<AtomPtr, double> AtomTarget;
typedef std::map<BondPtr, AtomTarget> BondEffects;
typedef std::map<BondPtr, double> BondCorrel;
typedef std::map<BondPtr, BondCorrel> BondBondCC;

class FlexGlobal;

typedef struct
{
	int index;
	int degree;
} BondDegree;

/**
 * \class FlexLocal
 * \brief Determines effect of each bond on all atom B factors, clusters them
 * and chooses some for refinement. */

class FlexLocal
{
public:
	FlexLocal();
	~FlexLocal();
	
	/** Sets polymer for which this FlexLocal object is responsible for,
	 * and the applicable magnitude of trial kick to worry about (not too
	 * much, not too small */
	void setPolymer(PolymerPtr pol, double shift);

	/** Evaluates target function of bond kicking parameters against data */
	static double getScore(void *object);
	
	/** Performs preliminary work and refinement simultaneously */
	void refine();
	
	/** Get current value of the magnitude of trial kick set */
	double getShift()
	{
		return _shift;
	}
	
	void setGetterSetter(Getter getter, Setter setter)
	{
		_getter = getter;
		_setter = setter;
	}
	
	void setWhacking(bool whack);
private:
	void createAtomTargets();
	AtomTarget currentAtomValues();
	void createClustering();
	void reorganiseBondOrder();
	double bondRelationship(BondPtr bi, BondPtr bj);
	void scanBondParams();
	void reflex();
	void clear();
	void propagateWhack();
	void setBondParam(BondPtr b, double k);
	double getBondParam(BondPtr b);
	double bondAtomCorrel(BondPtr b);
	double actualAtomChange(AtomPtr a);
	double targetForAtom(AtomPtr a);
	double getTotalBChange();
	
	double directSimilarity();
	static double sgetTotalBChange(void *object);
	double getTotalB();

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
	
	FlexGlobal *_flexGlobal;
	bool _useTarget;
	bool _usingWhack;
	int _afterBond;
	double _startB;
	double _threshold;
	double _increment;
	double _anchorB;
	double _negMult;
	int _direct;
	int _window;
	int _run;
	double _shift;
	Getter _getter;
	Setter _setter;

};

#endif


