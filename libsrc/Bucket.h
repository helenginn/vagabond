//
//  Bucket.h
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Bucket__
#define __vagabond__Bucket__

#include <stdio.h>
#include "shared_ptrs.h"
#include "MDNode.h"
#include <vector>
#include <map>
#include "vec3.h"
#include "../libccp4/csymlib.h"
#include "../libccp4/ccp4_spg.h"
#include "Plucker.h"

#define MAX_CHECK_DISTANCE 6.0
#define MIN_CHECK_DISTANCE 1.0

typedef std::map<double, Plucker *> PluckerMap;
typedef PluckerMap::iterator PluckerItr;

/**
 * \class Bucket
 * \brief Abstract class for adding solvent (via implemented method),
 * applying symmetry operations and scaling.
 */

class Bucket
{
public:
	Bucket();
	virtual void addSolvent() = 0;
	
	virtual ~Bucket()
	{
		for (PluckerItr it = _pluckerMap.begin(); it != _pluckerMap.end();
		     it++)
		{
			delete it->second;
			it->second = NULL;
		}
		
		delete _mdnode;
	}

	int getRandomValues(double left, double *right, double *angle);
	int getReallyRandomValues(double left, double *right, double *angle);

	void scaleSolvent();
	static double scaleSolventScore(void *object);
	double scaleAndAddSolventScore();
	void applySymOps(CSym::CCP4SPG *spaceGroup);
	void fourierTransform(int dir);
	void writeMillersToFile(std::string prefix, double maxRes);
	void abandonCalculations();
	
	void setCrystal(CrystalPtr crystal)
	{
		_crystal = crystal;
	}
	
	void setData(DiffractionPtr data)
	{
		_data = data;
	}
	
	static double getSolvScale(void *object)
	{
		return static_cast<Bucket *>(object)->_solvScale;
	}
	
	static void setSolvScale(void *object, double value)
	{
		static_cast<Bucket *>(object)->_solvScale = value;
	}

	static double getSolvBFac(void *object)
	{
		return static_cast<Bucket *>(object)->_solvBFac;
	}
	
	static void setSolvBFac(void *object, double value)
	{
		static_cast<Bucket *>(object)->_solvBFac = value;
	}

	void setupNodes(int split);
	void analyseSolvent(double distance);
	void loadAnalysis(std::string filename);
	
	bool isSolvent(vec3 pos);
	void processMaskedRegions();
	
	bool isSolvent(int index);
	Atom *nearbyAtom(int index);
	
	FFTPtr getMaskedRegions()
	{
		return _maskedRegions;
	}	
protected:
	void populateHistogram(MDNode *node, vec3 centre, vec3 left);
	int addAnalysisForSolventPos(MDNode *node, vec3 centre, double distance);
	CrystalWkr _crystal;
	MDNode *_mdnode;
	PluckerMap _pluckerMap;
	FFTPtr _solvent;
	FFTPtr _maskedRegions;
	DiffractionPtr _data;
	std::vector<Atom *> _atomPtrs;
	
	CrystalPtr getCrystal()
	{
		return _crystal.lock();
	}
private:
	/* Mask regions with protein = 0, solvent = 1 and protein/solvent
	* interface = 2 */
	double _solvScale;
	double _solvBFac;
	int _wanted;
	
	double _averages[3];
};

#endif /* defined(__vagabond__Bucket__) */
