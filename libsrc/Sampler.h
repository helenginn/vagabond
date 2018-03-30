//
//  Sampler.h
//  vagabond
//
//  Created by Helen Ginn on 25/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Sampler__
#define __vagabond__Sampler__

#include <stdio.h>
#include <vector>
#include "mat3x3.h"
#include "shared_ptrs.h"
#include "RefinementStrategy.h"
#include <string>
#include <map>

/* More of an abstraction, but will take a series of (bond) parameters,
* take a target function, and supply them to a refinement strategy. */

typedef enum
{
	RefinementBroad = 0,
	RefinementFine = 1,
	RefinementSidechain = 2,
	RefinementModelRMSDZero = 3,
	RefinementModelPos = 4,
	RefinementFlexibility = 5,
} RefinementType;

typedef enum
{
	ScoreTypeCorrel = 0,
	ScoreTypeMultiply = 1,
	ScoreTypeRFactor = 2,
	ScoreTypeModelRMSDZero = 3,
	ScoreTypeModelPos = 4,
	ScoreTypeModelFlexiness = 5,
} ScoreType;

typedef std::map<ParamOptionType, double> ParamMap;

class Sampler
{
public:
	Sampler();

	void addSampled(AtomPtr atom);
	void addSampled(std::vector<AtomPtr> atoms);

	void addTorsion(BondPtr bond, double range, double interval);
	void addTorsionBlur(BondPtr bond, double range, double interval);
	void addDampening(BondPtr bond, double range, double interval);
	void addBondLength(BondPtr bond, double range, double interval);
	void addBendAngle(BondPtr bond, double range, double interval);
	void addOccupancy(BondPtr bond, double range, double interval);

	/** Add sampled atoms from a given atom group and conformer name */
	void addSampledAtoms(AtomGroupPtr group, std::string conformer = "");
	void addRamachandranAngles(PolymerPtr polymer, int from, int to);
	void addAbsolutePosition(AbsolutePtr abs, double range, double interval);
	void addAbsoluteBFactor(AbsolutePtr abs, double range, double interval);
	void addRotamer(Sidechain *side, double range, double interval);
	void addMagicAngle(BondPtr bond, double range, double interval);
	void setCrystal(CrystalPtr crystal);
	bool sample(bool clear = true);

	void setupGrid();
	void setupNelderMead();

	void reportInDegrees()
	{
		_strategy->reportInDegrees();
	}

	void setCycles(int cycles)
	{
		_strategy->setCycles(cycles);
	}

	static double score(void *object)
	{
		return static_cast<Sampler *>(object)->getScore();
	}

	void setJobName(std::string newJob)
	{
		_jobName = newJob;
	}

	std::string getJobName()
	{
		return _jobName;
	}

	void setSilent(bool silent = true)
	{
		_strategy->setSilent(silent);
		_silent = silent;
	}

	void setVerbose(bool verbose = true)
	{
		_strategy->setVerbose(verbose);
	}

	int sampleSize()
	{
		return _sampled.size();
	}

	void setMock()
	{
		_mock = true;
	}

	void setScoreType(ScoreType type)
	{
		_scoreType = type;
	}

	void clearParams()
	{
		_params.clear();
	}

	void addParamType(ParamOptionType type, double value)
	{
		_params[type] = value; 
	}

	size_t paramCount()
	{
		return _params.size();
	}

	void copyParams(SamplerPtr sampler)
	{
		sampler->_params = _params;
	}


protected:
	BondPtr setupTorsionSet(BondPtr bond, int k, int bondNum,
	                        double range, double interval,
	bool addAngle = false, bool addFlex = false);
	BondPtr setupThoroughSet(BondPtr bond, int k, int bondNum,
	                         double range, double interval,
							 bool addAngle = false, bool addFlex = false);
	FFTPtr _fft;
	mat3x3 _real2Frac;

	int _refinedMagicAxisCount;
	virtual bool shouldRefineMagicAxis(BondPtr) { return false; }
	virtual double getScore();
private:
	void addAtomsForBond(BondPtr bond, int k);
	void addParamsForBond(BondPtr bond);
	void setupCloseAtoms();
	CrystalPtr _crystal;

	std::vector<AtomPtr> _sampled;
	std::vector<AtomPtr> _unsampled;
	std::vector<BondPtr> _bonds;
	ParamMap _params;

	bool _mock;
	bool _silent;

	std::string _jobName;
	ScoreType _scoreType;

	RefinementStrategyPtr _strategy;
};


#endif /* defined(__vagabond__Sampler__) */
