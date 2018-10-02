#ifndef __vagabond__MapScoreWorkspace__
#define __vagabond__MapScoreWorkspace__

#include "shared_ptrs.h"
#include "mat3x3.h"
#include <vector>

/**
 * \struct MapScoreWorkspace
 * \brief Stores temporary information for repeated calculations of
 * correlation between calculated and comparison density */

#define COORDVAL_FULL

typedef struct
{
#ifdef COORDVAL_FULL
	vec3 pos;
	double mask;
#endif
	double fo;
	double fc;
} CoordVal;

/** Scoring functions against map or PDB file */
typedef enum
{
	ScoreTypeCorrel = 0, /** Correlation between map and model density */
	ScoreTypeMultiply = 1, /** Weighted (by model) average of map voxels */
	ScoreTypeRFactor = 2, /** R factor in real space for electron density */
	ScoreTypeModelRMSDZero = 3, /** All ensemble positions against PDB */
	ScoreTypeBFactorAgreement = 4, /** All ensemble B factors against PDB */
	ScoreTypeModelPos = 5, /** Average ensemble position against PDB */
	ScoreTypeRMSDZero = 6, /** Sum of squares of anisotropic tensor */
	ScoreTypeAddDensity = 7, /** Sum of map voxels if model positive */
	ScoreTypeAddVoxels = 8, /** Number of map voxels if model positive */
	ScoreTypeHappiness = 9, /** Like CC but a reward function */
} ScoreType;

/* Don't forget - power of 2... */
typedef enum
{
	MapScoreFlagNone = 0,
	MapScoreFlagDifference = 1,
	MapScoreFlagReplaceWithObs = 2,
	MapScoreFlagNoSurround = 4,
	MapScoreFlagSkipScore = 8,
} MapScoreFlag;

typedef struct 
{
	ScoreType scoreType;
	CrystalPtr crystal;
	std::vector<AtomPtr> selectAtoms;
	std::vector<AtomPtr> extra;
	FFTPtr segment;
	FFTPtr fcSegment;
	FFTPtr constant;
	vec3 ave;
	mat3x3 basis;
	unsigned int flag;
} MapScoreWorkspace;


#endif
