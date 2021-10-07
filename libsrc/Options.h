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
#include "vagcout.h"
#include "Crystal.h"
#include "Notifiable.h"

typedef enum
{
	ModelFilePDB,
	ModelFileVagabond,
} ModelFile;

class StartScreen;
class Blob;

class Options
{
	friend class StartScreen;
public:
	Options(int argc, const char **argv);
	virtual ~Options()
	{

	}
	void run();
	void clear();
	void makeCout();

	static void setRuntimeOptions(OptionsPtr pointer)
	{
		Options::options = pointer;
		pointer->setManual(true);
	}

	static OptionsPtr getRuntimeOptions(bool make = true);

	void setNotify(Notifiable *notifiable)
	{
		_notify = notifiable;
		
		if (_filter != NULL)
		{
			_filter->setNotify(_notify);
		}
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
	
	void addCrystal(ParserPtr cryst)
	{
		CrystalPtr crystal = ToCrystalPtr(cryst);
		crystals.push_back(crystal);
	}

	static void loadBlobInGUI(Blob *b);
	
	static void changeSamplesAndFit(void *object, double n);

	static CrystalPtr getActiveCrystal()
	{
		if (!getRuntimeOptions()->crystals.size()) return CrystalPtr();
		return getRuntimeOptions()->crystals[_active];
	}
	
	static void setActiveCrystal(int c)
	{
		_active = c;
	}
	
	void removeLastCrystal()
	{
		crystals.pop_back();
		_active = 0;
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

	static double minRes()
	{
		return _minRes;
	}

	static double maxRes()
	{
		return _maxRes;
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
	
	static void setProbeRadius(void *, double r)
	{
		_probeRadius = r;
	}
	
	static double getShrink()
	{
		return _shrink;
	}
	
	static double getProbeRadius()
	{
		return _probeRadius;
	}
	
	static bool fitSolventParameters()
	{
		return _fitBucket;
	}
	
	static void setNSamples(void *, int n)
	{
		_nSamples = n;
	}

	static double getBMult()
	{
		return _bMult;
	}

	static double getBSubt()
	{
		return _bSubt;
	}
	
	static void setBSubt(void *, double bSubt)
	{
		_bSubt = bSubt;
		Options::getRuntimeOptions()->applyBMultiplier();
	}

	static void setBMult(void *, double bMult)
	{
		_bMult = bMult;
		Options::getRuntimeOptions()->applyBMultiplier();
	}

	static int getAddSolvent()
	{
		return _solvent;
	}
	
	/** Return which bond angles should be subject to refinement.
	 * Return values: 0 - no bond angle refinement
	 * 	              1 - Cb carbons 
	 * 	              2 - Cg angles for aromatics
	 * 	              3 - + glycine bond angles
	 * 	              4 - + various side chain bond angles
	 */
	static int getBondAngles()
	{
		return _bondAngles;
	}
	
	static bool getPeptideMovement()
	{
		return _peptideMovement;
	}
	
	static bool getRigidBody()
	{
		return _rRigidBody;
	}
	
	static bool isRefiningPositions()
	{
		return _rPosition;
	}

	static bool isRefiningInterflex()
	{
		return _rInter;
	}
	
	static bool isRefiningIntraflex()
	{
		return _rIntra;
	}
	
	static bool isAngling()
	{
		return _angling;
	}
	
	static void setSamplingFraction(void *object, double f)
	{
		_samplingFraction = f;	
	}
	
	static void setProteinSampling(void *object, double sampling)
	{
		_sampling = sampling;	
	}

	static double getSamplingFraction()
	{
		return _samplingFraction;
	}
	
	static double getProteinSampling()
	{
		return _sampling;
	}
	
	static int getMapType()
	{
		return _map;
	}
	
	static int getMaxRotations()
	{
		return _maxRot;
	}

	void setManual(bool manual)
	{
		_manual = manual;
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
	
	static double getUnmodelledFraction()
	{
		return _unmodelled;
	}
	
	static double getGlobalBFactor()
	{
		return _bReal;
	}
	
	static bool getScrew()
	{
		return _screw;
	}
	
	static bool getLowResMode()
	{
		return _lowResMode;
	}
	
	static std::string getLabSigF()
	{
		return _labSIGFP;
	}
	
	static std::string getLabF()
	{
		return _labFP;
	}
	
	static void setLabPhase(std::string phase)
	{
		_labPHI = phase;
	}
	
	static std::string getLabPhase()
	{
		return _labPHI;
	}

	static std::string getLabFWT()
	{
		return _labFWT;
	}
	
	static std::string getLabFree()
	{
		return _labFree;
	}
	
	static void setLabF(std::string f)
	{
		_labFP = f;
	}
	
	static void setLabFree(std::string f)
	{
		_labFree = f;
	}
	
	static void setLabSigF(std::string f)
	{
		_labSIGFP = f;
	}
	
	static void setGlobalBFactor(void *, double val)
	{
		_bReal = val;
	}
	
	static void resetGlobalBFactor()
	{
		_bReal = -1;
	}
	
	static int threads()
	{
		return _threads;
	}
	
	static void flagDensityChanged()
	{
		Options::getRuntimeOptions()->renderDensity();
	}
	
	static std::string anchorString()
	{
		return _anchor;
	}
	
	static bool usePartial()
	{
		return _usePartial;
	}
	
	static void pauseGUIFishing(bool on);
	virtual void focusOnPosition(vec3 pos, double dist = 25);
	void chelate();
	void renderWarp();
	void renderDensity();
	void disableDensityUpdate();
	static void statusMessage(std::string message, bool std_out = true);
	void agreementSummary();
	void applyBMultiplier();
	CrystalPtr openModel(std::string pdbName);
	void openMTZ(std::string mtzName);
	void recalculateFFT(bool saveState = true);
	void openInCoot();
	void findDisulphides();
	void adjustBFactor();
	void redirectToStream();

	static std::string rTypeString(RefinementType type);
private:
	void executeProtocol();
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
	bool _couted;
	int _globalCount;
	bool _tie;
	bool _manual;
	
	std::string _modelFile;
	std::string _mtzFile;
	std::string _addPDB;
	std::string _updatePDB;

	static bool _lowResMode;
	static bool _diagnostics;
	static bool _useRFree;
	static bool _usePartial;
	static bool _fitBucket;
	static ScalingType _scaleType;
	static std::string _solventFile;
	static double _kick;
	static int _solvent;
	static int _bondAngles;
	static int _nCycles;
	static int _pCycles;
	static int _threads;
	static bool _peptideMovement;
	static double _unmodelled;
	static double _bMult;
	static double _bSubt;
	static double _bReal;
	static double _bStart;
	static double _sampling;
	static double _samplingFraction;
	static int _nSamples;
	static int _map;
	static int _maxRot;
	static double _probeRadius;
	static double _shrink;
	std::string _diffMatrix;
	std::string _outputDir;
	static double _minRes;
	static double _maxRes;
	static std::string _anchor;
	static std::string _labFP;
	static std::string _labFWT;
	static std::string _labSIGFP;
	static std::string _labPHI;
	static std::string _labFree;
	
	static bool _rRigidBody;
	static bool _rPosition;
	static bool _rSidechains;
	static bool _refine;
	static bool _far;
	static bool _rInter;
	static bool _rIntra;
	static bool _convertWaters;
	static bool _hydrogens;
	static bool _screw;
	static bool _angling;
	static bool _refineGlobalParams;

	/*crystal*/
	static int _active;
	
	/* how many processes are currently locking GUI controls */
	int _processes;
	
	vagcout<char> *_filter;
};

#endif /* defined(__vagabond__Options__) */
