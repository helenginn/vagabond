//
//  WaterCluster.h
//  vagabond
//
//  Created by Helen Ginn on 01/07/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#ifndef __vagabond__WaterCluster__
#define __vagabond__WaterCluster__

#include <stdio.h>
#include "shared_ptrs.h"
#include "Molecule.h"
#include "Chromosomal.h"
#include <vector>
#include <map>
#include "Options.h"

typedef enum
{
	RestraintNone,
	RestraintRepel,
	RestraintHBond
} RestraintType;

typedef std::vector<RestraintType> Restraints;

typedef struct
{
	Restraints options;
	unsigned int picked;
} RestraintChoice;

typedef struct
{
	AtomPtr neighbour;
	bool close;
	RestraintChoice choice;
} RestraintMap;

typedef struct
{
	AtomPtr water;
	std::vector<RestraintMap> neighbours;
} WaterNeighboursPair;

/**
 * \class WaterCluster
 * \brief A group of atoms which includes water molecules and its
 * neighbouring atoms. Will perform hydrogen bonding duties.
 * 
 * At the moment this is mostly defunct.
 */

class WaterCluster : public Molecule, public Chromosomal
{
public:
	WaterCluster();
	WaterCluster(WaterCluster &other);
	virtual ~WaterCluster() {}
	
	virtual void addAtom(AtomPtr atom);

	virtual std::string getClassName()
	{
		return "WaterCluster";
	}
	
	void findNeighbours();
	
	static double score(void *object)
	{
		WaterCluster *cluster = static_cast<WaterCluster *>(object);
		return cluster->evaluateRestraints();
	}

	double scoreAgainstDensity();
	void refine();

	size_t waterCount()
	{
		return _waters.size();
	}
private:	
	int _modifySample;
	std::vector<AtomPtr> _waters;
	std::vector<WaterNeighboursPair> _pairs;

	double recalculateWaters();
	double evaluateRestraints();
	double evaluateRestraint(int sample, int i);
	void wipeBonding();
	void resetModels();
	
	/* Genetic things */
	virtual void randomise(double frac = 1.1);
	virtual double evaluate();
	virtual void mutate();
	virtual ChromosomalPtr makeCopy();
	virtual void geneticCode();
	virtual void haveSexWith(Chromosomal *_other);
	virtual void copyOver(Chromosomal *_other);
};


#endif


