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
#include "MapScoreWorkspace.h"
#include "mat3x3.h"
#include "shared_ptrs.h"
#include "RefinementStrategy.h"
#include <string>
#include <map>

/* More of an abstraction, but will take a series of (bond) parameters,
* take a target function, and supply them to a refinement strategy.
* Make sure the crystal target is set. */

/** Flags to set refinement strategies for a protein chain. */
typedef enum
{
	RefinementFine = 1, /** Refinement against electron density */
	RefinementSidechain = 2, /** Sidechains against electron density */
	RefinementModelRMSDZero = 3, /** Squeeze protein chain to PDB positions */
	RefinementModelPos = 4, /** Positions to PDB positions */
	RefinementRMSDZero = 5, /** Squeeze protein chain to minimise flexibility */
} RefinementType; 


typedef std::map<ParamOptionType, double> ParamMap;

class Sampler
{
public:
	Sampler();

	/**
	* 	 Add an atom from the sensitive area in real space
	*/
	void addSampled(AtomPtr atom);
	
	/**
	*  Add a bunch of atoms from the sensitive area in real space	
	*/
	void addSampled(std::vector<AtomPtr> atoms);

	/**
	* Add torsion angle to sampled parameters
	* \param bond bond to vary torsion around
	* \param range step size for angle in degrees
	* \param interval unused usually, but could be stopping point. 	 
	*/
	void addTorsion(BondPtr bond, double range, double interval);

	/**
	* Add torsion blur (kick function) to sampled parameters
	* \param bond bond to vary kick for
	* \param range step size in degree offset per degree
	* \param interval unused usually, but could be stopping point. 	 
	*/
	void addTorsionBlur(BondPtr bond, double range, double interval);
	void addDampening(BondPtr bond, double range, double interval);
	void addBondLength(BondPtr bond, double range, double interval);
	void addBendAngle(BondPtr bond, double range, double interval);
	void addOccupancy(BondPtr bond, double range, double interval);

	/** Add sampled atoms from a given atom group and conformer name */
	void addSampledAtoms(AtomGroupPtr group, std::string conformer = "");
	
	/** Add absolute position for an atom controlled by Absolute class */
	void addAbsolutePosition(AbsolutePtr abs, double range, double interval);
	/** Add absolute B factor for an atom controlled by Absolute class */
	void addAbsoluteBFactor(AbsolutePtr abs, double range, double interval);
	void addRotamer(Sidechain *side, double range, double interval);
	void addMagicAngle(BondPtr bond, double range, double interval);
	void setCrystal(CrystalPtr crystal);
	bool sample(bool clear = true);

	void setupGrid();
	void setupNelderMead();
	void setupStepSearch();

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

	/**
	* 	If mock is set, then the refinement is performed only to see the
	* 	change in target function and then the original values for the
	* 	parameters are restored.
	*/
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


	void addCustomParameter(void *object, Getter getter, Setter setter,
                                 double range, double interval,
                                 std::string name);
protected:
	BondPtr setupTorsionSet(BondPtr bond, int k, int bondNum,
	                        double range, double interval,
	bool addAngle = false, bool addFlex = false);
	BondPtr setupThoroughSet(BondPtr bond, int bondNum,
	                         double range, double interval,
							 bool addAngle = false, bool addFlex = false);
	FFTPtr _fft;
	mat3x3 _real2Frac;

	int _refinedMagicAxisCount;
	virtual bool shouldRefineMagicAxis(BondPtr) { return false; }
	virtual double getScore();
	void setupCloseAtoms();
	void setupScoreWithMap();
private:
	double preScanParameter(BondPtr bond, Getter getter, Setter setter,
	                      double stepSize);
	void addAtomsForBond(BondPtr bond, int k);
	void addParamsForBond(BondPtr bond);
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
	MapScoreWorkspace _workspace;
};


#endif /* defined(__vagabond__Sampler__) */
