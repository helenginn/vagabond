#ifndef __vagabond__MapScoreWorkspace__
#define __vagabond__MapScoreWorkspace__

#include "shared_ptrs.h"
#include "mat3x3.h"
#include <vector>

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
	ScoreTypeMultiply = 1, /** Weighted (by model) sum of map voxels */
	ScoreTypeRFactor = 2, /** R factor in real space for electron density */
	ScoreTypeModelRMSDZero = 3, /** All ensemble positions against PDB */
	ScoreTypeBFactorAgreement = 4, /** All ensemble B factors against PDB */
	ScoreTypeModelPos = 5, /** Average ensemble position against PDB */
	ScoreTypeRMSDZero = 6, /** Sum of squares of anisotropic tensor */
} ScoreType;

typedef enum
{
	MapScoreFlagNone = 0,
	MapScoreFlagDifference = 1,
	MapScoreFlagSubtractFc = 2,
} MapScoreFlag;

typedef struct 
{
	ScoreType scoreType;
	CrystalPtr crystal;
	std::vector<AtomPtr> selectAtoms;
	std::vector<AtomPtr> extra;
	FFTPtr segment;
	FFTPtr constant;
	vec3 ave;
	mat3x3 basis;
	unsigned int flag;
} MapScoreWorkspace;


#endif
