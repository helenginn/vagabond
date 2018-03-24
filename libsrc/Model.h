//
//  Model.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Model__
#define __vagabond__Model__

#include "shared_ptrs.h"
#include <stdio.h>
#include <vector>
#include "Distributor.h"
#include <mutex>
#include "Parser.h"

typedef struct
{
	mat3x3 basis;   /* Defines bond axis of previous bond */
	vec3 start;     /* position of last minor */
	vec3 old_start; /* position of torsion-defining atom */
	double torsion; /* Defines torsion of next atom */
	double occupancy;
} BondSample;

// Anything which is capable of predicting electron positions.
//

class Model : public boost::enable_shared_from_this<Model>, public Distributor, public Parser
{
public:
	Model();

	virtual FFTPtr getDistribution(bool quick = false, int n = -1) = 0;

	virtual void addToMonomer(MonomerPtr monomer);
	virtual void addToMolecule(MoleculePtr) {};

	virtual std::string getClassName() = 0;

	/* Static position if no blurring factors applied (for bonds) */

	/* Actual mean position of blurred positions (may not be same as static) */
	virtual vec3 getAbsolutePosition()
	{
		return _absolute;
	}

	/* Get blurred position array */
	virtual std::vector<BondSample> *getManyPositions() = 0;

	std::vector<vec3> polymerCorrectedPositions();
	virtual std::vector<BondSample> getFinalPositions();

	virtual double getEffectiveOccupancy() { return 1; }

	virtual double getMeanSquareDeviation() = 0;
	virtual mat3x3 getRealSpaceTensor();
	double biggestStdevDim();
	int fftGridLength();

	bool hasMolecule()
	{
		return !_molecule.expired();
	}

	MoleculePtr getMolecule()
	{
		return _molecule.lock();
	}

	void setMolecule(MoleculePtr mole)
	{
		_molecule = mole;
	}

	FFTPtr getZeroDistribution();
	virtual void propagateChange(int depth = -1, bool refresh = false);

	bool isBond()
	{
		return (getClassName() == "Bond");
	}

	bool isAbsolute()
	{
		return (getClassName() == "Absolute");
	}

	bool isAnchor()
	{
		return (getClassName() == "Anchor");
	}
	
	void setPolymerChanged()
	{
		_recalcFinal = true;	
	}

	vec3 longestAxis();
	std::vector<vec3> fishPositions();
protected:
	mat3x3 _realSpaceTensor;

	/* Molecule which can provide offsets/rotations/etc. */
	MoleculeWkr _molecule;

	/* What should be returned when asking for an atom's position
	* for drawing into a map... */
	vec3 _absolute;
	/* And a record of the final positions just for GUI display */
	std::vector<vec3> _finalPositions;

	/* Expect interference from GUI */
	static bool _useMutex;

	vec3 _longest;
	double _anisotropyExtent;

	virtual void getAnisotropy(bool withKabsch);
	double anisotropyExtent(bool withKabsch = false);
	double _isotropicAverage;
	bool _recalcFinal;
	std::vector<BondSample> _finalSamples;

	virtual std::string getParserIdentifier()
	{
		return "model"; 
	}

	virtual void addProperties();
	virtual void addObject(ParserPtr, std::string) {};
private:
	std::mutex guiLock;
};

#endif /* defined(__vagabond__Model__) */
