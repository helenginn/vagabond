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
#include "Crystal.h"

class Options
{
public:
	Options(int argc, const char **argv);

	static void setRuntimeOptions(OptionsPtr pointer)
	{
		Options::options = pointer;
	}

	static OptionsPtr getRuntimeOptions()
	{
		return options;
	}

	void run();

	size_t crystalCount()
	{
		return crystals.size();
	}

	CrystalPtr getCrystal(int i)
	{
		return crystals[i];
	}
	
private:
	static OptionsPtr options;

	void parse();
	void outputCrystalInfo();
	void refinementCycle(MoleculePtr molecule, int *count,
						 RefinementType type);
	void refineAll(RefinementType type, int numCycles, int *count);

	std::vector<std::string> arguments;

	std::vector<ObjectPtr> objects;
	std::vector<CrystalPtr> crystals;
	std::vector<DatasetPtr> datasets;
	std::vector<DiffractionPtr> diffractions;

	int _numCycles;
	bool _tie;
	std::string _outputDir;
};

#endif /* defined(__vagabond__Options__) */
