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

/* More of an abstraction, but will take a series of (bond) parameters,
 * take a target function, and supply them to a refinement strategy. */

typedef enum
{
	RefinementBroad = 0,
	RefinementFine = 1,
	RefinementFineBlur = 2,
	RefinementModelRMSD = 3,
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

class Sampler
{
public:
	Sampler();

	void addSampled(AtomPtr atom);
	void addSampled(std::vector<AtomPtr> atoms);

	void addOverallKickAndDampen(PolymerPtr polymer);
	void addSidechainDampen(PolymerPtr polymer);
	
	void addTorsion(BondPtr bond, double range, double interval);
	void addTorsionBlur(BondPtr bond, double range, double interval);
	void addDampening(BondPtr bond, double range, double interval);
	void addBondLength(BondPtr bond, double range, double interval);
	void addBendAngle(BondPtr bond, double range, double interval);
	void addOccupancy(BondPtr bond, double range, double interval);
	void addSampledBackbone(PolymerPtr polymer, int from = 0, int to = 0);
	void addSampledSidechains(PolymerPtr polymer);
	void addSampledAtoms(AtomGroupPtr group, std::string conformer = "");
	void addRamachandranAngles(PolymerPtr polymer, int from, int to);
	void addAbsolutePosition(AbsolutePtr abs, double range, double interval);
	void addAbsoluteBFactor(AbsolutePtr abs, double range, double interval);
	void addRotamer(Sidechain *side, double range, double interval);
//	void addMagicAngle(BondPtr bond, double range, double interval);
	void setCrystal(CrystalPtr crystal);
	bool sample(bool clear = true);

	void setupGrid();
	void setupNelderMead();
	void setupSnake();

	void setJointSampling()
	{
		_joint = true;
	}

	void copyTarget(Sampler *other)
	{
		_fft = other->_fft;
		_real2Frac = other->_real2Frac;
		_scoreType = other->_scoreType;
	}

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

	void setTargetFlexibility(double value)
	{
		_overallFlex = value;
	}

	std::vector<double> getNextResult(int num);

protected:
	BondPtr setupTorsionSet(BondPtr bond, int k, int bondNum,
							double range, double interval,
						 bool addDampening = false);
	FFTPtr _fft;
	mat3x3 _real2Frac;

	int _refinedMagicAxisCount;
	virtual bool shouldRefineMagicAxis(BondPtr bond) { return false; }
	virtual double getScore();
private:
	void addAtomsForBond(BondPtr bond, int k);

	std::vector<AtomPtr> _sampled;
	std::vector<AtomPtr> _unsampled;
	std::vector<BondPtr> _bonds;
	bool _mock;
	bool _joint;
	bool _silent;
	double _overallFlex;

	std::string _jobName;
	ScoreType _scoreType;

	RefinementStrategyPtr _strategy;
};


#endif /* defined(__vagabond__Sampler__) */
