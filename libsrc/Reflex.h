//
//  Reflex.h
//  vagabond
//
//  Created by Helen Ginn, 2018
//  Copyright © 2018 Helen Ginn. All rights reserved.
//

#ifndef Reflex_h
#define Reflex_h

/**
 * \class Reflex
 * \brief Assigns some kind of flex estimation to regions of a protein.
 * Great Flexpectations! */

#include "shared_ptrs.h"
#include "MapScoreWorkspace.h"
#include "RefineMat3x3.h"

typedef enum
{
	PieceTypeUnassigned,
	PieceTypeAtom,
	PieceTypeResidue,
} PieceType;

class Reflex
{
public:
	Reflex();
	
	/** Perform the calculations to assign B factor estimations */
	void calculate();

	/** Set the polymer for B factor estimations */
	void setPolymer(PolymerPtr pol)
	{
		_polymer = pol;
	}
	
	/** Set number of "pieces" (e.g. residues or atoms) to be considered in a
	 * single unit when calculating the local flexstimation */
	void setPieceCount(size_t count)
	{
		_pieceCount = count;
	}
	
	static double getBFactor(void *object)
	{
		return static_cast<Reflex *>(object)->_bFactor;
	}
	
	static void setBFactor(void *object, double value)
	{
		static_cast<Reflex *>(object)->_bFactor = value;
	}

	static double getAbsoluteScale(void *object)
	{
		return static_cast<Reflex *>(object)->_kAbs;
	}
	
	static void setAbsoluteScale(void *object, double value)
	{
		static_cast<Reflex *>(object)->_kAbs = value;
	}

private:
	void prepareSegments(AtomGroupPtr segment, bool preprocess = true);
	void segmentPolymer();
	void refineFlex();
	static double bFactorScore(void *object);
	
	size_t _pieceCount;
	PieceType _pieceType;
	PolymerPtr _polymer;
	MapScoreWorkspace _workspace;
	
	FFTPtr _obs, _calc;
	int _round;
	double _bFactor;
	double _kAbs;
	double _intercept;
	double _gradient;
	double _calcSize;
	RefineMat3x3Ptr _aniso;
	
	double _calcScaleFactor;
	
	std::vector<AtomGroupPtr> _segments;
};

#endif
