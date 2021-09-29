//
//  FlexLocal.h
//  vagabond
//
//  Created by Helen Ginn, 2018
//  Copyright © 2018 Helen Ginn. All rights reserved.
//

#ifndef FlexLocal_h
#define FlexLocal_h

class FlexLocal;
#include "shared_ptrs.h"
#include <hcsrc/RefinementStrategy.h>
#include "MapScoreWorkspace.h"
#include <map>

typedef std::map<AtomPtr, double> AtomTarget;
typedef std::map<BondPtr, AtomTarget> BondEffects;
typedef std::map<BondPtr, double> BondCorrel;
typedef std::map<BondPtr, BondCorrel> BondBondCC;

#include "SVDBond.h"

class FlexGlobal;

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
	
	/** Refines chain-dependent multipliers */
	void refineChainMults(AnchorPtr anch);

	/** Get current value of the magnitude of trial kick set */
	double getShift()
	{
		return _shift;
	}
	
	bool didChange()
	{
		return _changed;
	}
	
	PolymerPtr getPolymer()
	{
		return _polymer;
	}
	
	SVDBond *getSVD()
	{
		return _svd;
	}
private:
	void svd();
	void findAtomsAndBonds();
	void refineClusters();
	void scanBondParams();
	void clear();
	void propagateWhack();
	void setBondParam(BondPtr b, double w, double k);
	
	PolymerPtr _polymer;
	AtomGroupPtr _bb;

	BondEffects _bondEffects;
	
	std::vector<AtomPtr> _atoms;
	std::vector<BondPtr> _bonds;
	SVDBond *_svd;
	
	FlexGlobal *_flexGlobal;
	int _run;
	double _shift;

	bool _prepared;
	bool _changed;
	MapScoreWorkspace _workspace;
};

#endif


