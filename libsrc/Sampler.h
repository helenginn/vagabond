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
	RefinementFine = 1
} RefinementType;

class Sampler
{
public:
	Sampler();

	void addSampled(AtomPtr atom);
	void addTorsion(BondPtr bond, double range, double interval);
	void addTorsionBlur(BondPtr bond, double range, double interval);
	void addTorsionNextBlur(BondPtr bond, double range, double interval);
	void addBondLength(BondPtr bond, double range, double interval);
	void addBendBlur(BondPtr bond, double range, double interval);
	void addBendAngle(BondPtr bond, double range, double interval);
	void addOccupancy(BondPtr bond, double range, double interval);
	void addAbsolutePosition(AbsolutePtr abs, double range, double interval);
	void addAbsoluteBFactor(AbsolutePtr abs, double range, double interval);
	void setCrystal(CrystalPtr crystal);
	void sample();

	void setupGrid();
	void setupNelderMead();

	void setJointSampling()
	{
		_joint = true;
	}

	void reportInDegrees()
	{
		_strategy->reportInDegrees();
	}

	static double score(void *object)
	{
		return static_cast<Sampler *>(object)->getScore();
	}

	void setJobName(std::string newJob)
	{
		_jobName = newJob;
	}

	int sampleSize()
	{
		return _sampled.size();
	}

	void setMock()
	{
		_mock = true;
	}
protected:
	void setupDoubleTorsion(BondPtr bond, int k, int bondNum, int resNum,
							double range, double interval);

private:
	double getScore();

	std::vector<AtomPtr> _sampled;
	std::vector<AtomPtr> _unsampled;
	std::vector<BondPtr> _bonds;
	FFTPtr _fft;
	bool _mock;
	mat3x3 _real2hkl;
	bool _joint;

	std::string _jobName;

	RefinementStrategyPtr _strategy;
};


#endif /* defined(__vagabond__Sampler__) */
