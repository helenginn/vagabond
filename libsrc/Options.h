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

	static double getKick()
	{
		return _kick;
	}

	static double getDampen()
	{
		return _dampen;
	}

	static int enableTests()
	{
		return _enableTests;
	}

	static double getBStart()
	{
		return _bStart;
	}
private:
	static OptionsPtr options;

	void parse();
	void displayHelp();
	void outputCrystalInfo();
	void refinementCycle(MoleculePtr molecule, int *count,
						 RefinementType type);
	void refineAll(RefinementType type, int numCycles, int *count,
				   bool keepGoing = true);
	bool parseJoke(std::string arg);

	std::vector<std::string> arguments;

	std::vector<ObjectPtr> objects;
	std::vector<CrystalPtr> crystals;
	std::vector<DatasetPtr> datasets;
	std::vector<DiffractionPtr> diffractions;

	int _numCycles;
	bool _tie;
	static double _kick;
	static double _dampen;
	static int _enableTests;
	static double _bStart;
	std::string _outputDir;
};

#endif /* defined(__vagabond__Options__) */
