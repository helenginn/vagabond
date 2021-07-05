#ifndef __vagabond__MapScoreWorkspace__
#define __vagabond__MapScoreWorkspace__

#include "shared_ptrs.h"
#include <hcsrc/mat3x3.h>
#include <hcsrc/Timer.h>
#include <vector>
#include <map>

/**
 * \struct MapScoreWorkspace
 * \brief Stores temporary information for repeated calculations of
 * correlation between calculated and comparison density */

//#define COORDVAL_FULL

typedef struct
{
#ifdef COORDVAL_FULL
	vec3 pos;
#endif
	double fo;
	double fc;
	double weight;
} CoordVal;

/** Scoring functions against map or PDB file */
typedef enum
{
	ScoreTypeCorrel = 0, /** Correlation between map and model density */
	ScoreTypeMultiply = 1, /** Weighted (by model) average of map voxels */
	ScoreTypeRFactor = 2, /** R factor in real space for electron density */
	ScoreTypeModelPos = 5, /** Average ensemble position against PDB */
	ScoreTypeAddDensity = 7, /** Sum of map voxels if model positive */
	ScoreTypeAddVoxels = 8, /** Number of map voxels if model positive */
	ScoreTypeHappiness = 9, /** Like CC but a reward function */
	ScoreTypeZero = 10, /** Will always return 0 */
	ScoreTypeCentroid = 11, /** Difference of centroid against PDB */
	ScoreTypeCopyToSmaller = 12, /** Difference of centroid against PDB */
	ScoreTypeMouse = 13, /** Difference to Mouse position */
	ScoreTypeSavedPos = 14, /** To previously saved atom position */
} ScoreType;

/* Don't forget - power of 2... */
typedef enum
{
	MapScoreFlagNone = 0,
	MapScoreFlagDifference = 1,
	MapScoreFlagReplaceWithObs = 2,
	MapScoreFlagNoSurround = 4,
	MapScoreFlagSkipScore = 8,
	MapScoreFlagPosOnly = 16,
	MapScoreFlagNegOnly = 32,
} MapScoreFlag;

typedef enum
{
	MapScoreTypeNone,
	MapScoreTypeCorrel,
	MapScoreTypeCorrelCopy,
	MapScoreTypeRadialMagnitude,
	MapScoreTypeCopyToSmaller,
	MapScoreAddNoWrap,
} MapScoreType;

class MapScoreWorkspace
{
public:
	ScoreType scoreType;
	CrystalPtr crystal;
	AtomGroupPtr selectAtoms;
	AtomGroupPtr extra;
	VagFFTPtr comparison;
	VagFFTPtr segment;
	VagFFTPtr recip;
	std::vector<CoordVal> vals;
	std::string filename;
	unsigned int flag;
	bool recalc;
	Timer *tBonds;
	Timer *tMap;
	Timer *tScore;
};

inline void setup_space(MapScoreWorkspace *w)
{
	w->scoreType = ScoreTypeCorrel;
	w->crystal = CrystalPtr();
	w->selectAtoms = AtomGroupPtr();
	w->extra = AtomGroupPtr();
	w->comparison = VagFFTPtr();
	w->segment = VagFFTPtr();
	w->recip = VagFFTPtr();
	w->vals.clear();
	w->filename = "";
	w->flag = MapScoreFlagNone;
	w->recalc = false;
	w->tBonds = new Timer("updating bonds", false);
	w->tMap = new Timer("populating map", false);
	w->tScore = new Timer("scoring function", false);
}

#endif
