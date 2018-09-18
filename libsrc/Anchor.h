//
//  Anchor.h
//  vagabond
//
//  Created by Helen Ginn 2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#ifndef __vagabond__Anchor__
#define __vagabond__Anchor__

#include "ExplicitModel.h"

/**
 * \class Anchor
 * \brief Model for an anchor point in a tied-up polymer. */

class Anchor : public ExplicitModel
{
public:
	Anchor(AbsolutePtr absolute);
	Anchor();
	virtual ~Anchor() {};

	std::vector<BondSample> *getManyPositions(void *object = NULL);

	/** Returns the offsets for an anchor residue on which a molecule may
	* calculate translations and offsets.
	* \return x, y, z values centred around the origin.
	* */
	std::vector<vec3> getSphereAngles()
	{
		return _sphereAngles;
	}

	void setNeighbouringAtoms(AtomPtr nPre, AtomPtr nAtom, 
	                          AtomPtr cAtom, AtomPtr cPost);

	AtomPtr getOtherAtom(AtomPtr calling);

	void setOccupancies(std::vector<double> occ)
	{
		_occupancies = occ;
	}
	
	double getBFactor()
	{
		return _bFactor;
	}
	
	void setBFactor(double b)
	{
		_bFactor = b;
	}
	
	AtomPtr getAtom()
	{
		return _atom.lock();
	}

	static void setTransTensor11(void *object, double value)
	{
		static_cast<Anchor *>(object)->_translation.vals[0] = value;
		static_cast<Anchor *>(object)->propagateChange(-1, true);
	}

	static void setTransTensor12(void *object, double value)
	{
		static_cast<Anchor *>(object)->_translation.vals[1] = value;
		static_cast<Anchor *>(object)->_translation.vals[3] = value;
		static_cast<Anchor *>(object)->propagateChange(-1, true);
	}

	static void setTransTensor21(void *object, double value)
	{
		static_cast<Anchor *>(object)->_translation.vals[3] = value;
		static_cast<Anchor *>(object)->propagateChange(-1, true);
	}

	static void setTransTensor13(void *object, double value)
	{
		static_cast<Anchor *>(object)->_translation.vals[2] = value;
		static_cast<Anchor *>(object)->_translation.vals[6] = value;
		static_cast<Anchor *>(object)->propagateChange(-1, true);
	}

	static void setTransTensor31(void *object, double value)
	{
		static_cast<Anchor *>(object)->_translation.vals[6] = value;
		static_cast<Anchor *>(object)->propagateChange(-1, true);
	}

	static void setTransTensor22(void *object, double value)
	{
		static_cast<Anchor *>(object)->_translation.vals[4] = value;
		static_cast<Anchor *>(object)->propagateChange(-1, true);
	}

	static void setTransTensor32(void *object, double value)
	{
		static_cast<Anchor *>(object)->_translation.vals[7] = value;
		static_cast<Anchor *>(object)->propagateChange(-1, true);
	}

	static void setTransTensor23(void *object, double value)
	{
		static_cast<Anchor *>(object)->_translation.vals[5] = value;
		static_cast<Anchor *>(object)->_translation.vals[7] = value;
		static_cast<Anchor *>(object)->propagateChange(-1, true);
	}

	static void setTransTensor33(void *object, double value)
	{
		static_cast<Anchor *>(object)->_translation.vals[8] = value;
		static_cast<Anchor *>(object)->propagateChange(-1, true);
	}

	static double getTransTensor11(void *object)
	{
		return static_cast<Anchor *>(object)->_translation.vals[0];
	}

	static double getTransTensor21(void *object)
	{
		return static_cast<Anchor *>(object)->_translation.vals[3];
	}

	static double getTransTensor12(void *object)
	{
		return static_cast<Anchor *>(object)->_translation.vals[1];
	}

	static double getTransTensor31(void *object)
	{
		return static_cast<Anchor *>(object)->_translation.vals[6];
	}

	static double getTransTensor13(void *object)
	{
		return static_cast<Anchor *>(object)->_translation.vals[2];
	}

	static double getTransTensor22(void *object)
	{
		return static_cast<Anchor *>(object)->_translation.vals[4];
	}

	static double getTransTensor23(void *object)
	{
		return static_cast<Anchor *>(object)->_translation.vals[5];
	}

	static double getTransTensor32(void *object)
	{
		return static_cast<Anchor *>(object)->_translation.vals[7];
	}

	static double getTransTensor33(void *object)
	{
		return static_cast<Anchor *>(object)->_translation.vals[8];
	}
	
	static double getRotCentreX(void *object)
	{
		return static_cast<Anchor *>(object)->_rotCentre.x;
	}
	
	static double getRotCentreY(void *object)
	{
		return static_cast<Anchor *>(object)->_rotCentre.y;
	}
	
	static double getRotCentreZ(void *object)
	{
		return static_cast<Anchor *>(object)->_rotCentre.z;
	}
	
	static void setRotCentreX(void *object, double value)
	{
		static_cast<Anchor *>(object)->_rotCentre.x = value;
		static_cast<Anchor *>(object)->propagateChange(-1, true);
	}
	
	static void setRotCentreY(void *object, double value)
	{
		static_cast<Anchor *>(object)->_rotCentre.y = value;
		static_cast<Anchor *>(object)->propagateChange(-1, true);
	}
	
	static void setRotCentreZ(void *object, double value)
	{
		static_cast<Anchor *>(object)->_rotCentre.z = value;
		static_cast<Anchor *>(object)->propagateChange(-1, true);
	}

	static double getRotVecX(void *object)
	{
		return static_cast<Anchor *>(object)->_rotVec.x;
	}
	
	static double getRotVecY(void *object)
	{
		return static_cast<Anchor *>(object)->_rotVec.y;
	}
	
	static double getRotVecZ(void *object)
	{
		return static_cast<Anchor *>(object)->_rotVec.z;
	}
	
	static void setRotVecX(void *object, double value)
	{
		static_cast<Anchor *>(object)->_rotVec.x = value;
		static_cast<Anchor *>(object)->propagateChange(-1, true);
	}
	
	static void setRotVecY(void *object, double value)
	{
		static_cast<Anchor *>(object)->_rotVec.y = value;
		static_cast<Anchor *>(object)->propagateChange(-1, true);
	}
	
	static void setRotVecZ(void *object, double value)
	{
		static_cast<Anchor *>(object)->_rotVec.z = value;
		static_cast<Anchor *>(object)->propagateChange(-1, true);
	}

	virtual void propagateChange(int depth = -1, bool refresh = false);
	virtual std::string shortDesc();
protected:
	virtual std::string getClassName()
	{
		return "Anchor";
	}
	
	virtual std::string getParserIdentifier();

	virtual void sanityCheck();
	virtual void addProperties();
	virtual void linkReference(ParserPtr object, std::string category); 

private:
	void createStartPositions(Atom *callAtom);
	void translateStartPositions();
	void rotateBases();
	mat3x3 _translation;
	vec3 _rotVec;
	vec3 _rotCentre;
	vec3 _position;

	double _bFactor;
	AtomWkr _atom;
	AtomWkr _nAtom, _cAtom, _nPre, _cPost;
	
	vec3 _nDir, _nDir2;
	vec3 _cDir, _cDir2;

	std::vector<vec3> _sphereAngles;
	std::vector<double> _occupancies;
};

#endif
