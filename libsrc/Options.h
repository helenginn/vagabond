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

class StartScreen;

class Options
{
	friend class StartScreen;
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
	int parseParameter(std::string arg, std::string prefix, std::string *ptr);
	int parseParameter(std::string arg, std::string prefix, bool *ptr);
	int parseParameter(std::string arg, std::string prefix, int *ptr);
	
	static bool ignoreRFree()
	{
		return !_useRFree;
	}
	
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

	static double maxRes()
	{
		return _maxRes;
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

	static int getNSamples()
	{
		return _nSamples;
	}
	
	template <class T>
	static void setNSamples(T n)
	{
		_nSamples = n;
	}

	static double getBMult()
	{
		return _bMult;
	}

	static void setBMult(double bMult)
	{
		_bMult = bMult;
		Options::getRuntimeOptions()->applyBMultiplier();
	}

	static int getAddSolvent()
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

	static ScalingType getScalingType()
	{
		return _scaleType;
	}
	
	static void setScalingType(ScalingType shell)
	{
		_scaleType = shell;
	}
	
	static bool makeDiagnostics()
	{
		return _diagnostics;
	}
	
	static double getGlobalBFactor()
	{
		return _bReal;
	}
	
	static void setGlobalBFactor(double val)
	{
		_bReal = val;
	}
	
	static void resetGlobalBFactor()
	{
		_bReal = -1;
	}
	
	static void flagDensityChanged()
	{
		Options::getRuntimeOptions()->renderDensity();
	}
	
	static std::string anchorString()
	{
		return _anchor;
	}
	
	void omitScan();
	void reflex();
	void scanBondParams();
	void renderDensity();
	void statusMessage(std::string message);
	void agreementSummary();
	void previousState();
	void refineAll(RefinementType type, int numCycles, int *count = NULL,
	               bool keepGoing = false);
	void applyBMultiplier();
	void openModel(std::string pdbName);
	void openMTZ(std::string mtzName);
	void recalculateFFT(bool saveState = true);
	void openInCoot();
	void fitWholeMolecule(bool translation, bool rotation);
	void findDisulphides();
	void adjustBFactor();

	static std::string rTypeString(RefinementType type);
private:
	void executeScript();
	static OptionsPtr options;
	Notifiable *_notify;
	void notifyGUI(bool enable);

	void writeCommandLine();
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

	bool _parsed;
	int _globalCount;
	bool _tie;
	bool _manual;
	std::string _scriptName;
	
	std::string _modelFile;
	std::string _mtzFile;

	static bool _diagnostics;
	static bool _useRFree;
	static ScalingType _scaleType;
	static bool _powder;
	static std::string _solventFile;
	static double _kick;
	static int _solvent;
	static double _dampen;
	static double _bMult;
	static double _bReal;
	static int _enableTests;
	static double _bStart;
	static double _sampling;
	static int _nSamples;
	std::string _diffMatrix;
	std::string _outputDir;
	static double _minRes;
	static double _maxRes;
	static std::string _anchor;
};

#endif /* defined(__vagabond__Options__) */
