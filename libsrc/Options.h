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
#include "Notifiable.h"

typedef enum
{
	ModelFilePDB,
	ModelFileVagabond,
} ModelFile;


class Options
{
public:
	Options(int argc, const char **argv);
	void run();

	double getActiveCrystalDStar();
	
	static void setRuntimeOptions(OptionsPtr pointer)
	{
		Options::options = pointer;
		pointer->setManual(true);
	}

	static OptionsPtr getRuntimeOptions()
	{
		return options;
	}

	void setNotify(Notifiable *notifiable)
	{
		_notify = notifiable;
	}

	Notifiable *getNotify()
	{
		return _notify;
	}

	size_t crystalCount()
	{
		return crystals.size();
	}

	CrystalPtr getCrystal(int i)
	{
		return crystals[i];
	}

	CrystalPtr getActiveCrystal()
	{
		if (!crystals.size()) return CrystalPtr();
		return crystals[0];
	}

	DiffractionPtr getActiveData()
	{
		if (!diffractions.size()) return DiffractionPtr();
		return diffractions[0];
	}

	int parseParameter(std::string arg, std::string prefix, double *ptr);
	
	static double getKick()
	{
		return _kick;
	}

	static double getDampen()
	{
		return _dampen;
	}

	static double minRes()
	{
		return _minRes;
	}

	static int enableTests()
	{
		return _enableTests;
	}

	static double getBStart()
	{
		return _bStart;
	}

	static std::string getSolventFile()
	{
		return _solventFile;
	}


	static double getBMult()
	{
		return _bMult;
	}

	static void setBMult(double bMult)
	{
		_bMult = bMult;
	}

	static bool getAddSolvent()
	{
		return _solvent;
	}
	
	static void setProteinSampling(double sampling)
	{
		_sampling = sampling;	
	}

	static double getProteinSampling()
	{
		return _sampling;
	}

	void setManual(bool manual)
	{
		_manual = manual;
	}

	static bool shouldPowder()
	{
		return _powder;
	}

	void statusMessage(std::string message);
	void agreementSummary();
	void previousState();
	void backboneAnalysis();
	void refineAll(RefinementType type, int numCycles, int *count = NULL,
	               bool keepGoing = false);
	void superimposeAll(CrystalPtr crystal = CrystalPtr());
	void applyBMultiplier();
	void openModel(std::string pdbName);
	void openMTZ(std::string mtzName);
	void recalculateFFT(bool saveState = true);
	void fitWholeMolecule(bool translation, bool rotation);
	void findDisulphides();

	static std::string rTypeString(RefinementType type);
private:
	void executeScript();
	static OptionsPtr options;
	Notifiable *_notify;
	void notifyGUI(bool enable);

	void parse();
	void displayHelp();
	void outputCrystalInfo();
	void refinementCycle(MoleculePtr molecule, RefinementType type);
	bool parseJoke(std::string arg);
	void diffMatrix();

	std::vector<std::string> arguments;

	std::vector<ObjectPtr> objects;
	std::vector<CrystalPtr> crystals;
	std::vector<DatasetPtr> datasets;
	std::vector<DiffractionPtr> diffractions;

	int _numCycles;
	int _globalCount;
	bool _tie;
	bool _manual;
	std::string _scriptName;

	static bool _powder;
	static std::string _solventFile;
	static double _kick;
	static int _solvent;
	static double _dampen;
	static double _bMult;
	static int _enableTests;
	static double _bStart;
	static double _sampling;
	std::string _diffMatrix;
	std::string _outputDir;
	static double _minRes;
};

#endif /* defined(__vagabond__Options__) */
