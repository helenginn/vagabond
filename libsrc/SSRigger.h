//
//  SSRigger.h
//  vagabond
//
//  Created by Helen Ginn on 18/07/2017.
//  Copyright (c) 2018 Strubi. All rights reserved.
//

#include "shared_ptrs.h"

class SSRigger
{
public:
	SSRigger();	
	
	void setCrystal(CrystalPtr crystal)
	{
		_crystal = crystal;
	}
	
	void findDisulphides();
private:
	void findCysteineSulphurs();
	void findCloseCysteines();
	void convertCysteine(AtomPtr oneAtom);
	
	CrystalPtr _crystal;	
	std::vector<AtomPtr> _cysSGs;
};
