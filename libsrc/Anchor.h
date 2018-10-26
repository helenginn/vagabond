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
	
	void addWhack(WhackPtr whack)
	{
		_whacks.push_back(whack);
	}
	
	size_t whackCount()
	{
		return _whacks.size();
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
protected:
	virtual std::string getClassName()
	{
		return "Anchor";
	}
	
	virtual std::string getParserIdentifier();
	virtual void createStartPositions(Atom *callAtom);

	virtual void sanityCheck();
	virtual void addProperties();
	virtual void linkReference(ParserPtr object, std::string category); 
	virtual void addObject(ParserPtr object, std::string category);

	double _bFactor;
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
	
	vec3 _nDir, _nDir2;
	vec3 _cDir, _cDir2;

	std::map<Atom *, SamplePair> _samples;
	std::vector<vec3> _sphereAngles;
	std::vector<double> _occupancies;
	std::vector<WhackPtr> _whacks;
	
	bool _disableWhacks;
};

#endif
