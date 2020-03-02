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
#include "Quat4Refine.h"

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

	/** Mean position of blurred positions without including whole-molecule
	* 	deviations. Should not be used for finding final atom position
	* 	calculations. */
	virtual std::vector<BondSample> *getManyPositions(void *object = NULL) = 0;

	/** Positions and associated data including whole-molecule deviations.
	* 	Will return from cache if not flagged to recalculate. */
	virtual std::vector<BondSample> getFinalPositions();
	
	vec3 *samplePointer(int i)
	{
		return &_storedSamples[i].start;
	}
	
	unsigned int sampleSize()
	{
		return _storedSamples.size();
	}

	vec3 getSpecificPosition(int i)
	{
		if (i > (int)_finalPositions.size()) return empty_vec3();
		return getFinalPositions()[i].start;
	}
	
	void writePositionsToFile(std::string filename);
	
	void refreshPositions();

	std::vector<vec3> fishPositions(vec3 *average = NULL);

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
	
	void clearTwists();
	
	size_t positionCount()
	{
		return _storedSamples.size();
	}

	void addTwist(TwistPtr twist)
	{
		_twists.push_back(twist);
	}
	
	size_t twistCount()
	{
		return _twists.size();
	}
	
	void addShifts(RefinementStrategyPtr strategy,
                              AtomGroupPtr clearGroup);
	
	void setRotCentre(vec3 centre)
	{
		_centre = centre;
	}

	virtual bool canFish()
	{
		return true;
	}
	
	double getLowestZ()
	{
		return _lowestZ;
	}

	static double propagate(void *obj);
	static std::vector<vec3> makeCloud(double totalPoints, 
	                                   double b,
	                                   std::vector<double> &occs);
protected:
	bool _changedSamples;
	virtual void sanityCheck();

	std::vector<TwistPtr> _twists;

	virtual void getAnisotropy(bool withKabsch = false);
	virtual double anisotropyExtent(bool withKabsch = false);
	std::vector<BondSample> _storedSamples;
	vec3 meanOfManyPositions(std::vector<BondSample> *positions);

	/** Will define torsion basis as:
	* x: along line of 0º torsion angle.
	* y: completes the right-handed coordinate system
	* z: along bond direction, from heavy-to-light alignment atoms.
	*/
	mat3x3 makeTorsionBasis(vec3 hPos, vec3 maPos,
	                        vec3 miPos, vec3 lPos, double *newAngle = NULL);
	std::mutex guiLock;
	
	/* lowest Z value for re-ordering */
	double _lowestZ;
private:
	std::vector<BondSample> _finalSamples;
	int _modifySample;
	vec3 _centre;
	
	Quat4Refine _shift;
	Quat4Refine _rotation;
};

#endif
