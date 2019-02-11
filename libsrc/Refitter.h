//
//  Reflex.h
//  vagabond
//
//  Created by Helen Ginn, 2018
//  Copyright © 2018 Helen Ginn. All rights reserved.
//

#ifndef Refitter_h
#define Refitter_h

#include "shared_ptrs.h"
#include "Sampler.h"

typedef struct
{
	ParamBandPtr list;
	BondPtr last;
	double score;
} BondTrial;

class Refitter : public Sampler
{
public:

	/** Constructor for refitter
	 * @param bond parent-level bond after split
	 * @param forwards indicates on C-terminal side of anchor. */

	Refitter(BondPtr bond, BondPtr target, bool forwards);
	
	void refit();

private:
	void makePeptides180();
	BondPtr scanBond(ParamBandPtr trial, BondPtr former, 
                           std::vector<BondTrial> *populate);
	BondPtr generateCandidates(ParamBandPtr trial, int residue,
	                           std::vector<BondTrial> *populate);
	BondPtr findNextBond(BondPtr bond);

	void pruneUnreasonable(double distance = -1, bool landing = false);
	void strictPrune(double angle = 30);
	void removeDuplicates();
	void addSolutions();
	void applyTrial(ParamBandPtr trial, bool propagate = true,
	                int advance = 0);
	double accumulativeBondLength(BondPtr current);
	void refineCentroid(BondTrial *trial, int tries = 10);
	void collectTorsionsIntoTrial(BondTrial *trial);
	
	void removeLanding();
	
	ParamBandPtr collectTorsionsIntoSolution();
	BondPtr advanceBondConformer(BondPtr bond, int advance);
	
	std::vector<BondTrial> _trialList;
	std::vector<ParamBandPtr> _solutions;

	AtomGroupPtr _group;
	BondPtr _bond;
	BondPtr _target;
	PolymerPtr _polymer;
	
	/* forwards direction towards C terminus */
	bool _forwards;
	
	int _mStart, _mEnd, _mLand;
	
	double _degstep;

	/* atom name for other part of peptide bond */
	std::string _pepAtom;
};

#endif


