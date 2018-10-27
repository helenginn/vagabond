//
//  ExplicitModel.h
//  vagabond
//
//  Created by Helen Ginn 2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#ifndef __vagabond__ExplicitModel__
#define __vagabond__ExplicitModel__

#include "Model.h"

/**
 * \class ExplicitModel
* \brief Abstract class providing a template for models which rely on
* explicit atom position calculations.
*
* These will include the Bond models and the anchor position for a polymer.
 */

/** \struct BondSample
 *  \brief Transfers information between Bonds for calculation of positional
 *  information.
 */
typedef struct
{
	mat3x3 basis;     /**< Defines bond axis of previous bond */
	vec3 start;       /**< position of last minor */
	vec3 old_start;   /**< position of torsion-defining atom */
	double torsion;   /**< To be filled in by next bond temporarily */
	double kickValue; /**< Contains number between 0 and 1, kick multiplier */
	double occupancy; /**< Relative occupancy (usually 1) */
} BondSample;

class ExplicitModel : public Model
{
public:
	ExplicitModel();
	virtual ~ExplicitModel() {};

	virtual bool hasExplicitPositions()
	{
		return true;
	}
	FFTPtr makeDistribution();

	/** Mean position of blurred positions without including whole-molecule
	* 	deviations. Should not be used for finding final atom position
	* 	calculations. */
	virtual std::vector<BondSample> *getManyPositions(void *object = NULL) = 0;

	/** Positions and associated data including whole-molecule deviations.
	* 	Will return from cache if not flagged to recalculate. */
	virtual std::vector<BondSample> getFinalPositions();

	vec3 getSpecificPosition(int i)
	{
		if (i > _finalPositions.size()) return empty_vec3();
		return getFinalPositions()[i].start;
	}
	
	void writePositionsToFile(std::string filename);
	
	void refreshPositions();

	std::vector<vec3> fishPositions(vec3 *average = NULL);

	void setPosN(int choice, double value);
	double getPosN(int choice);

	static void setPosX(void *object, double x)
	{
		ExplicitModel *expl = static_cast<ExplicitModel *>(object);
		expl->setPosN(0, x);
	}

	static void setPosY(void *object, double y)
	{
		ExplicitModel *expl = static_cast<ExplicitModel *>(object);
		expl->setPosN(1, y);
	}

	static void setPosZ(void *object, double z)
	{
		ExplicitModel *expl = static_cast<ExplicitModel *>(object);
		expl->setPosN(2, z);
	}

	static double getPosZ(void *object)
	{
		ExplicitModel *expl = static_cast<ExplicitModel *>(object);
		return expl->getPosN(2);
	}

	static double getPosX(void *object)
	{
		ExplicitModel *expl = static_cast<ExplicitModel *>(object);
		return expl->getPosN(0);
	}

	static double getPosY(void *object)
	{
		ExplicitModel *expl = static_cast<ExplicitModel *>(object);
		return expl->getPosN(1);
	}

	void setModifiedSample(int i)
	{
		_modifySample = i;
	}
	
	void clearModifiedSample()
	{
		_modifySample = -1;
	}


	/** Returns the B factor (function is a misnomer). */
	virtual double getMeanSquareDeviation();

	virtual vec3 longestAxis();
	
	static void useMutex()
	{
		_useMutex = true;
	}

	virtual size_t downstreamBondCount(int)
	{
		return 1;
	}
	
	virtual int downstreamBondNum(Bond *, int *)
	{
		return 0;
	}
	
	bool changedSamples()
	{
		return _changedSamples;
	}
protected:
	FFTPtr makeRealSpaceDistribution();
	void addRealSpacePositions(FFTPtr real, vec3 offset);
	bool _changedSamples;
	virtual void sanityCheck();

	virtual void getAnisotropy(bool withKabsch = false);
	virtual double anisotropyExtent(bool withKabsch = false);
	std::vector<BondSample> _storedSamples;

	/** Will define torsion basis as:
	* x: along line of 0º torsion angle.
	* y: completes the right-handed coordinate system
	* z: along bond direction, from heavy-to-light alignment atoms.
	*/
	mat3x3 makeTorsionBasis(vec3 hPos, vec3 maPos,
	                        vec3 miPos, vec3 lPos, double *newAngle = NULL);
private:
	std::mutex guiLock;
	std::vector<BondSample> _finalSamples;
	int _modifySample;
};

#endif
