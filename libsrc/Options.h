//
//  Options.h
//  vagabond
//
//  Created by Helen Ginn on 22/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Options__
#define __vagabond__Options__

#include <stdio.h>
#include <vector>
#include <string>
#include "shared_ptrs.h"

class Options
{
public:
	Options(int argc, const char **argv);

	void run();
	
private:
	void parse();
	void outputCrystalInfo();

	std::vector<std::string> arguments;

	std::vector<ObjectPtr> objects;
	std::vector<CrystalPtr> crystals;
	std::vector<DatasetPtr> datasets;
	std::vector<DiffractionPtr> diffractions;

};

#endif /* defined(__vagabond__Options__) */
