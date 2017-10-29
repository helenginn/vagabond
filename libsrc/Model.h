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

typedef struct
{
	mat3x3 basis;
	vec3 start;     /* position of last minor */
	vec3 old_start; /* position of torsion-defining atom */
	double torsion;
	double occupancy;
} BondSample;

typedef enum
{
	BondSampleThorough,
	BondSampleStatic,
	BondSampleMonteCarlo
} BondSampleStyle;

// Anything which is capable of predicting electron positions.
//

class Model : public std::enable_shared_from_this<Model>, public Distributor
{
public:
	Model();

	virtual FFTPtr getDistribution(bool quick = false) = 0;

	virtual void addToMonomer(MonomerPtr monomer);
	virtual void addToMolecule(MoleculePtr molecule);

	virtual std::string getClassName() = 0;

	/* Static position if no blurring factors applied (for bonds) */
	virtual vec3 getStaticPosition() = 0;

	/* Actual mean position of blurred positions (may not be same as static) */
	virtual vec3 getAbsolutePosition() = 0;

	/* Get blurred position array */
	virtual std::vector<BondSample> *getManyPositions(BondSampleStyle style) = 0;

	virtual std::vector<BondSample> getFinalPositions();

	virtual double getMeanSquareDeviation() = 0;
	virtual mat3x3 getRealSpaceTensor()
	{
		return _realSpaceTensor;
	}

	FFTPtr getZeroDistribution();
	virtual void propagateChange();

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
protected:
	mat3x3 _realSpaceTensor;
private:
	std::vector<AtomPtr> atoms;
};

#endif /* defined(__vagabond__Model__) */
