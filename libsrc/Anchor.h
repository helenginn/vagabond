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

class RefineMat3x3;

typedef struct
{
	std::vector<BondSample> samples;
	bool changed;
} SamplePair;

class Anchor : public ExplicitModel
{
public:
	/* Substantiate an Anchor using a pre-existing Absolute atom. Absolute
	 * atom ought to be a backbone nitrogen atom. */
	Anchor(AbsolutePtr absolute);
	
	/* Substantiate a default Anchor, to be used by Parser */
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

	/** To be called with the appropriate neighbouring atoms along
	 * the backbone.
	 * \param nPre two atoms prior to nitrogen (previous Calpha atom)
	 * \param nAtom one atom prior to nitrogen (previous carbon atom)
	 * \param cAtom one atom further on than nitrogen (following Calpha atom)
	 * \param cPost two atoms further than nitrogen (following carbon atom)
	 */
	void setNeighbouringAtoms(AtomPtr nPre, AtomPtr nAtom, 
	                          AtomPtr cAtom, AtomPtr cPost);

	/** Returns the opposing atom on the other side of the bond.
	 * \param calling Either the following C-alpha or preceding carbon atom.
	 * \return empty AtomPtr() if calling atom is not one of the options, or
	 * the opposite option if valid. */
	AtomPtr getOtherAtom(AtomPtr calling);
	
	AtomPtr getNAtom()
	{
		return _nAtom.lock();
	}
	
	AtomPtr getCAtom()
	{
		return _cAtom.lock();
	}

	void setOccupancies(std::vector<double> occ)
	{
		_occupancies = occ;
	}
	
	double getBFactor()
	{
		return _bFactor;
	}
	
	static double sgetBFactor(void *object)
	{
		Anchor *anch = static_cast<Anchor *>(object);
		return anch->_bFactor;
	}
	
	static void ssetBFactor(void *object, double b)
	{
		Anchor *anch = static_cast<Anchor *>(object);
		anch->_bFactor = b;
		anch->propagateChange(-1, true);
	}
	
	void setBFactor(double b)
	{
		_bFactor = b;
	}
	
	AtomPtr getAtom()
	{
		return _atom.lock();
	}

	void addWhack(WhackPtr whack)
	{
		_whacks.push_back(whack);
	}
	
	size_t whackCount()
	{
		return _whacks.size();
	}
	
	WhackPtr getWhack(int i)
	{
		return _whacks[i];
	}
	
	void addTranslationParameters(RefinementStrategyPtr strategy,
	                              double mult = 1);
	void addLibrationParameters(RefinementStrategyPtr strategy,
	                              double mult = 1);

	static void cleanup(void *object)
	{
		static_cast<Anchor *>(object)->propagateChange(-1, true);
	}

	void recalculateWhacks();
	virtual void propagateChange(int depth = -1, bool refresh = false);
	virtual std::string shortDesc();


	static double getPosX(void *object)
	{
		Anchor *anch = static_cast<Anchor *>(object);
		return anch->_position.x;
	}

	static double getPosY(void *object)
	{
		Anchor *anch = static_cast<Anchor *>(object);
		return anch->_position.y;
	}

	static double getPosZ(void *object)
	{
		Anchor *anch = static_cast<Anchor *>(object);
		return anch->_position.z;
	}

	static void setPosX(void *object, double x)
	{
		Anchor *anch = static_cast<Anchor *>(object);
		anch->_position.x = x;
		anch->propagateChange(-1);
	}

	static void setPosY(void *object, double y)
	{
		Anchor *anch = static_cast<Anchor *>(object);
		anch->_position.y = y;
		anch->propagateChange(-1);
	}

	static void setPosZ(void *object, double z)
	{
		Anchor *anch = static_cast<Anchor *>(object);
		anch->_position.z = z;
		anch->propagateChange(-1);
	}

	static double getAlpha(void *object)
	{
		Anchor *anch = static_cast<Anchor *>(object);
		return anch->_alpha;
	}

	static double getBeta(void *object)
	{
		Anchor *anch = static_cast<Anchor *>(object);
		return anch->_beta;
	}

	static double getGamma(void *object)
	{
		Anchor *anch = static_cast<Anchor *>(object);
		return anch->_gamma;
	}

	static void setAlpha(void *object, double val)
	{
		Anchor *anch = static_cast<Anchor *>(object);
		anch->_alpha = val;
		anch->propagateChange(-1);
	}

	static void setBeta(void *object, double val)
	{
		Anchor *anch = static_cast<Anchor *>(object);
		anch->_beta = val;
		anch->propagateChange(-1);
	}

	static void setGamma(void *object, double val)
	{
		Anchor *anch = static_cast<Anchor *>(object);
		anch->_gamma = val;
		anch->propagateChange(-1);
	}
	
	/* Set matrix describing principle axes of polymer. This should
	 * help choose sensible starting parameters for libration */
	void setPolymerBasis(mat3x3 basis);

	virtual std::string getClassName()
	{
		return "Anchor";
	}
protected:
	
	virtual std::string getParserIdentifier();
	virtual void createStartPositions(Atom *callAtom);

	virtual void sanityCheck();
	virtual void addProperties();
	virtual void linkReference(ParserPtr object, std::string category); 
	virtual void addObject(ParserPtr object, std::string category);

	double _bFactor;
	
	/* Euler angles for modifying _nAtom and _cAtom */
	double _alpha, _beta, _gamma;
	RefineMat3x3Ptr _trans;
	RefineMat3x3Ptr _libration;
	
	bool _changedN;
	bool _changedC;
	
	mat3x3 _rotation;
	vec3 _position;
	AtomWkr _atom;
	AtomWkr _nAtom, _cAtom;
private:
	void translateStartPositions();
	void rotateBases();
	void fixCentroid();
	mat3x3 getAnchorRotation();
	mat3x3 _libMotion;
	
	vec3 _nDir, _nDir2;
	vec3 _cDir, _cDir2;

	std::map<Atom *, SamplePair> _samples;
	std::vector<vec3> _sphereAngles;
	std::vector<double> _occupancies;
	std::vector<WhackPtr> _whacks;
	
	bool _disableWhacks;
};

#endif
