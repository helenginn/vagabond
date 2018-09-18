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
typedef std::map<AtomPtr, AtomTarget> AtomTargets;

class FlexLocal
{
public:
	FlexLocal();
	
	void setPolymer(PolymerPtr pol)
	{
		_polymer = pol;
	}

	void scanBondParams();
private:
	void createAtomTargets();

	PolymerPtr _polymer;
	AtomTarget _atomTargets;
	AtomTargets _pairwiseTargets;
	
	std::vector<AtomPtr> _atoms;
	std::vector<BondPtr> _bonds;

};

#endif


