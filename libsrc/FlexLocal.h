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
#include "RefinementStrategy.h"
#include "MapScoreWorkspace.h"
#include <map>

typedef std::map<AtomPtr, double> AtomTarget;
typedef std::map<BondPtr, AtomTarget> BondEffects;
typedef std::map<BondPtr, double> BondCorrel;
typedef std::map<BondPtr, BondCorrel> BondBondCC;

#include "SVDBond.h"

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

class SVDBond;

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
	
	/** Go button. Performs preliminary work and refinement */
	void refine();
	
	/** Get current value of the magnitude of trial kick set */
	double getShift()
	{
		return _shift;
	}
	
	/** Sets what the flexibility parameter is focusing on for refinement -
	 * default is to use kicks, but usually updated to using Whacks */
	void setGetterSetter(Getter getter, Setter setter)
	{
		_getter = getter;
		_setter = setter;
	}
private:
	void createAtomTargets();
	void createClustering();
	void svd();
	double bondRelationship(BondPtr bi, BondPtr bj);
	void scanBondParams();
	void reflex();
	void clear();
	void propagateWhack();
	void setBondParam(BondPtr b, double w, double k);
	double actualAtomChange(AtomPtr a);
	double targetForAtom(AtomPtr a);
	
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
	SVDBond *_svd;
	
	FlexGlobal *_flexGlobal;
	bool _useTarget;
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

	bool _prepared;
	MapScoreWorkspace _workspace;
};

#endif


