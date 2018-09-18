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
	PolymerPtr _polymer;
};


#endif


