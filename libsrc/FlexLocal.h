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

typedef struct
{
	AtomPtr atom;
	double value;
} AtomValuePair;

typedef std::map<AtomPtr, double> AtomTarget;
typedef std::map<BondPtr, AtomValuePair> BondEffects;

class FlexLocal
{
public:
	FlexLocal();
	
	void setPolymer(PolymerPtr pol)
	{
		_polymer = pol;
	}

	void scanBondParams();
	static double getScore(void *object);
	void refine();
private:
	void createAtomTargets();
	AtomTarget currentAtomValues();
	
	double directSimilarity();

	PolymerPtr _polymer;

	AtomTarget _atomTargets;
	AtomTarget _atomOriginal;
	BondEffects _bondEffects;
	
	/* Targets, not actual */
	std::vector<AtomPtr> _atoms;
	std::vector<BondPtr> _bonds;
	
	double _anchorB;
	int _direct;
	int _window;
	int _run;
	double _shift;

};

#endif


