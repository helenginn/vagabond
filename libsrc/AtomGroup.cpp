//
//  AtomGroup.cpp
//  vagabond
//
//  Created by Helen Ginn on 25/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#define MAP_VALUE_CUTOFF 20.

#include "Timer.h"
#include "AtomGroup.h"
#include <climits>
#include "Atom.h"
#include "Bond.h"
#include <sstream>
#include "Crystal.h"
#include <iomanip>
#include "CSV.h"
#include "maths.h"
#include "Plucker.h"
#include "Shouter.h"
#include "../libccp4/ccp4_spg.h"
#include "Options.h"
#include "fftw3d.h"
#include <time.h>

AtomPtr AtomGroup::findAtom(std::string atomType)
{
	for (size_t i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getAtomName() == atomType)
		{
			return atom(i);
		}
	}

	return AtomPtr();
}

AtomList AtomGroup::findAtoms(std::string atomType, int resNum)
{
	AtomList atoms = findAtoms(atomType);

	for (size_t i = 0; i < atoms.size(); i++)
	{
		if (atoms[i]->getResidueNum() != resNum)
		{
			atoms.erase(atoms.begin() + i);
			i--;
		}
	}

	return atoms;
}

AtomPtr AtomGroup::findAtom(std::string atomType, std::string confID)
{
	AtomList atoms = findAtoms(atomType);

	for (size_t i = 0; i < atoms.size(); i++)
	{
		if (atoms[i]->getAlternativeConformer() == confID)
		{
			return atoms[i];
		}
	}

	return AtomPtr();
}

void AtomGroup::addAtomsFrom(AtomGroupPtr group)
{
	for (size_t i = 0; i < group->atomCount(); i++)
	{
		addAtom(group->atom(i));	
	}
}

std::map<std::string, size_t> AtomGroup::conformerMap()
{
	std::map<std::string, size_t> conformerList;

	for (size_t i = 0; i < atomCount(); i++)
	{
		std::string conformer = atom(i)->getAlternativeConformer();

		if (!conformerList.count(conformer))
		{
			conformerList[conformer] = 0;
		}

		conformerList[conformer]++;
	}

	return conformerList;
}

int AtomGroup::conformerCount()
{
	std::map<std::string, size_t> conformerList = conformerMap();

	return conformerList.size();
}

std::string AtomGroup::conformer(size_t i)
{
	if (i > conformerMap().size()) return "";

	std::map<std::string, size_t> conformerList = conformerMap();
	std::map<std::string, size_t>::iterator it = conformerList.begin();

	for (size_t j = 0; j < i; j++) it++;

	return it->first;
}

AtomList AtomGroup::findAtoms(std::string atomType)
{
	AtomList list;

	for (size_t i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getAtomName() == atomType)
		{
			list.push_back(atom(i));
		}
	}

	return list;
}

double AtomGroup::totalElectrons()
{
	double total = 0;

	for (size_t i = 0; i < atomCount(); i++)
	{
		total += atom(i)->getElectronCount();
	}

	return total;
}

std::string AtomGroup::getPDBContribution(PDBType pdbType, CrystalPtr crystal,
                                          int conformer)
{
	std::ostringstream stream;
	size_t numConf = 0;

	if (!atomCount())
	{
		return "";
	}

	if (pdbType == PDBTypeEnsemble)
	{
		/* Give up if not explicit */
		if (!atom(0)->getModel()->hasExplicitPositions())
		{
			stream << atom(0)->averagePDBContribution(false, false);
			return stream.str();
		}

		/* Get the total number of conformers to worry about */
		std::vector<BondSample> *samples;
		samples = atom(0)->getExplicitModel()->getManyPositions();

		numConf = samples->size();

		for (size_t i = 0; i < atomCount(); i++)
		{
			if (atom(i)->getWeighting() <= 0)
			{
				continue;
			}

			stream << atom(i)->getPDBContribution(conformer);
		}

		return stream.str();
	}

	for (size_t i = 0; i < atomCount(); i++)
	{
		bool samePos = (pdbType == PDBTypeSamePosition);
		bool sameB = (pdbType == PDBTypeSameBFactor);
		stream << atom(i)->averagePDBContribution(samePos, sameB);

		if (crystal)
		{
			stream << atom(i)->anisouPDBLine(crystal);
		}
	}

	return stream.str();
}

double AtomGroup::getAverageDisplacement()
{
	double sum = 0;
	double count = 0;

	for (size_t i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getElectronCount() <= 1)
		{
			continue;
		}

		double val = atom(i)->posDisplacement();

		sum += val;
		count++;
	}

	return sum / count;
}

double AtomGroup::getAverageBFactor(bool initial)
{
	double sum = 0;
	double count = 0;

	for (size_t i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getElectronCount() <= 1)
		{
			continue;
		}

		if (initial)
		{
			sum += atom(i)->getInitialBFactor();
		}
		else
		{
			double val = atom(i)->getBFactor();
			sum += val;
		}

		count++;
	}

	return sum / count;
}

AtomGroup::AtomGroup()
{
	_beenTied = false;
	_timesRefined = 0;
	_largestNum = -INT_MAX;
}

void AtomGroup::propagateChange()
{
	for (size_t i = 0; i < atomCount(); i++)
	{
		atom(i)->getModel()->propagateChange(0);
	}
}

void AtomGroup::refreshPositions(bool quick)
{
	for (size_t i = 0; i < atomCount(); i++)
	{
		if (!atom(i)) continue;

		atom(i)->getModel()->propagateChange(0);
		atom(i)->getModel()->refreshPositions();
	}

	if (quick) return;

	AtomList list = topLevelAtoms();

	for (size_t i = 0; i < list.size(); i++)
	{
		AtomPtr atom = list[i].lock();
		atom->getModel()->propagateChange(-1, true);
	}
}


int AtomGroup::totalElectrons(int *fcWeighted)
{
	double sum = 0;
	double weighted = 0;

	for (size_t i = 0; i < atomCount(); i++)
	{
		double e = atom(i)->getElectronCount();
		sum += e;
		double weight = atom(i)->getWeighting();
		weighted += e * weight;
	}

	*fcWeighted = weighted;

	return sum;
}

void AtomGroup::setWeighting(double value)
{
	for (size_t i = 0; i < atomCount(); i++)
	{
		atom(i)->setWeighting(value);
	}
}

AtomList AtomGroup::topLevelAtoms()
{
	if (!atomCount()) return AtomList();

	AtomList list;
	
	for (size_t i = 0; i < conformerCount(); i++)
	{
		std::string conf = conformer(i);
		size_t j = 0;
		AtomPtr topAtom = atom(0);

		while (topAtom->getAlternativeConformer() != conf)
		{
			if (j >= atomCount())
			{
				goto giveup;
			}

			j++;
			topAtom = atom(j);
		}

		while (true)
		{
			if (!topAtom->getModel()->isBond())
			{
				break;
			}

			BondPtr bond = ToBondPtr(topAtom->getModel());

			if (!hasAtom(bond->getMajor()))
			{
				break;
			}

			if (shouldRefineAtom(bond->getMajor()))
			{
				topAtom = bond->getMajor();
			}
		}

		list.push_back(topAtom);

		giveup:
		continue;
	}

	return list;
}

bool AtomGroup::hasAtom(AtomPtr anAtom)
{
	bool found = false;

	if (!anAtom) return false;

	for (size_t i = 0; i < atomCount(); i++)
	{
		if (atom(i) == anAtom)
		{
			found = true;
		}
	}

	return found;
}

void AtomGroup::setTargetRefinement(CrystalPtr target, RefinementType rType)
{
	_target = target;
	_rType = rType;
}

void AtomGroup::privateRefine()
{
	time_t wall_start;
	time(&wall_start);
	std::cout << "Refining one residue." << std::endl;
	refine(_target, _rType);
	shout_timer(wall_start, "refinement");
}

void AtomGroup::addAtom(AtomPtr atom)
{
	std::vector<AtomPtr>::iterator it;
	it = std::find(_atoms.begin(), _atoms.end(), atom);

	if (it == _atoms.end())
	{
		_atoms.push_back(atom);
		if (atom->getAtomNum() > _largestNum)
		{
			_largestNum = atom->getAtomNum();
		}
	}
}

Plucker *AtomGroup::makePluckableWaters()
{
	Plucker *plucker = new Plucker();
	plucker->setGranularity(0.2);

	for (int i = 0; i < atomCount(); i++)
	{
		AtomPtr atm = atom(i);

		if (!atm->isHeteroAtom() || !(atm->getAtomName() == "O"))
		{
			continue;
		}

		// we have a water
		atm->cacheCloseWaters(4.);

		if (atm->pluckCount())
		{
			double occupancy = atm->getModel()->getEffectiveOccupancy();
			plucker->addPluckable(&*atm, occupancy);
		}
	}

	return plucker;
}


AtomPtr AtomGroup::getClosestAtom(CrystalPtr crystal, vec3 pos)
{
	double small_dist = FLT_MAX;
	AtomPtr best;

	for (int i = 0; i < atomCount(); i++)
	{
		vec3 tmp = atom(i)->getPositionInAsu();
		bool closeish = vec3_near_vec3_box(tmp, pos, small_dist);

		if (!closeish)
		{
			continue;
		}
		else
		{
			vec3 diff = vec3_subtract_vec3(tmp, pos);
			double dist = vec3_sqlength(diff);
			if (dist < small_dist)
			{
				small_dist = dist;
				best = atom(i);
			}
		}
	}

	return best;
}

void AtomGroup::refine(CrystalPtr target, RefinementType rType)
{
	AtomList topAtoms = topLevelAtoms();
	bool refineAngles = shouldRefineAngles();
	_timesRefined++;

	ScoreType scoreType = ScoreTypeModelPos;
	int maxTries = 0;

	switch (rType) {
		case RefinementModelPos:
		scoreType = ScoreTypeModelPos;
		maxTries = 60;
		break;

		case RefinementFine:
		scoreType = ScoreTypeCorrel;
		maxTries = 2;
		break;

		case RefinementModelRMSDZero:
		scoreType = ScoreTypeModelRMSDZero;
		maxTries = 10;
		break;

		case RefinementRMSDZero:
		scoreType = ScoreTypeRMSDZero;
		maxTries = 10;
		break;

		default:
		shout_at_helen("Unimplemented refinement option?");
		break;
	}
	
	while (topAtoms.size() > 0)
	{
		BondPtr topBond;
		int count = 0;

		bool changed = true;
		bool addFlex = (rType == RefinementFine);

		while (changed && count < maxTries)
		{
			setupNelderMead();
			setCrystal(target);
			setCycles(16);

			for (int i = 0; i < topAtoms.size(); i++)
			{
				if (topAtoms[i].expired())
				{
					continue;
				}
				
				AtomPtr topAtom = topAtoms[i].lock();
				
				if (!topAtom->getModel()->isBond())
				{
					continue;
				}
				
				BondPtr bond = ToBondPtr(topAtom->getModel());

				if (!bond->isRefinable())
				{
					continue;
				}

				if (i == 0)
				{
					setJobName("torsion_" +  bond->shortDesc());
				}

				if (rType != RefinementFine)
				{
					topBond = setupThoroughSet(bond, false);
				}
				else
				{
					topBond = setupThoroughSet(bond);
				}
			}

			setScoreType(scoreType);

			for (size_t l = 0; l < _includeForRefine.size(); l++)
			{
				addSampledAtoms(_includeForRefine[l]);
			}

			if (rType == RefinementModelPos 
			    || rType == RefinementFine 
			    || rType == RefinementModelRMSDZero
			    || rType == RefinementRMSDZero)
			{
				setSilent();
			}

			changed = sample();
			count++;
		}

		if (!topBond)
		{
			break;
		}

		AtomPtr topAtom = topBond->getMinor();
		topAtoms = findAtoms(topAtom->getAtomName());
		
		if (!hasAtom(topAtom))
		{
			break;
		}
	}

	_includeForRefine.clear();
	refreshPositions(true);
}

double AtomGroup::scoreWithMap(ScoreType scoreType, CrystalPtr crystal, bool plot, unsigned int flags)
{
	OptionsPtr options = Options::getRuntimeOptions();
	DiffractionPtr data = options->getActiveData();
	if (!data)
	{
		return 0;
	}

	std::vector<AtomPtr> selected;
	for (size_t i = 0; i < atomCount(); i++)
	{
		selected.push_back(atom(i));
	}
	
	MapScoreWorkspace workspace;
	workspace.scoreType = scoreType;
	workspace.crystal = crystal;
	workspace.selectAtoms = selected;
	workspace.segment = FFTPtr();
	workspace.ave = empty_vec3();
	workspace.basis = make_mat3x3();
	workspace.flag = flags;

	return scoreWithMapGeneral(&workspace, plot);
}

FFTPtr AtomGroup::prepareMapSegment(CrystalPtr crystal,
                                    std::vector<AtomPtr> selected,
                                    mat3x3 *basis, vec3 *ave)
{
	double maxDistance = 0;
	FFTPtr map = crystal->getFFT();
	mat3x3 real2Frac = crystal->getReal2Frac();

	vec3 sum = make_vec3(0, 0, 0);

	/* Find centroid of atom set */
	for (size_t i = 0; i < selected.size(); i++)
	{
		selected[i]->getModel()->refreshPositions();
		vec3 offset = selected[i]->getModel()->getAbsolutePosition();
		sum = vec3_add_vec3(sum, offset);	
	}

	vec3_mult(&sum, 1 / (double)selected.size());

	if (sum.x != sum.x)
	{
		return FFTPtr();
	}

	sum = crystal->snapToGrid(sum);
	*ave = sum;

	/* Find the longest distance from the centroid to determine
	* the max FFT dimensions.*/
	double ns[3];
	ns[0] = 0; ns[1] = 0; ns[2] = 0;

	for (size_t i = 0; i < selected.size(); i++)
	{
		/* Refresh absolute position */
		selected[i]->getModel()->refreshPositions();
		vec3 offset = selected[i]->getModel()->getAbsolutePosition();

		vec3 diff = vec3_subtract_vec3(offset, sum);

		ns[0] = std::max(fabs(diff.x), ns[0]);
		ns[1] = std::max(fabs(diff.y), ns[1]);
		ns[2] = std::max(fabs(diff.z), ns[2]);
	}
	
	double maxDStar = Options::getRuntimeOptions()->getActiveCrystalDStar();
	double scales = 1.0 / (2 * maxDStar);

	/** Adding a buffer region */
	long nl[3];
	const double buffer = 1.5;

	for (int i = 0; i < 3; i++)
	{
		nl[i] = 2 * (ns[i] + buffer);
	}
	
	/* Correction of non-orthogonal unit cells, I think */

	mat3x3 crystal_basis = crystal->getFFT()->getBasisInverse();
	mat3x3 angs = make_mat3x3();
	mat3x3_scale(&angs, nl[0], nl[1], nl[2]);
	mat3x3 resized = mat3x3_mult_mat3x3(crystal_basis, angs);

	for (int i = 0; i < 9; i+=3)
	{
		nl[i / 3] = std::max(std::max(resized.vals[i], resized.vals[i + 1]), 
		                     resized.vals[i + 2]);
	}

	/* Calculate appropriate box size and setup FFT */

	FFTPtr segment = FFTPtr(new FFT());
	segment->create(nl[0], nl[1], nl[2]);
	
	/* Basis of the segment itself needs to be in voxels to angstroms */
	*basis = crystal->getFFT()->getBasis();
	segment->setBasis(*basis);

	/* but the basis of the map workspace should be a mini-
	 * reciprocal unit cell */
	double toReal[3];

	for (int i = 0; i < 3; i++)
	{
		toReal[i] = 1 / (segment->scales[i] * nl[i]);
	}

	mat3x3_scale(basis, nl[0], nl[1], nl[2]);
	*basis = mat3x3_inverse(*basis);

	return segment;
}

double AtomGroup::scoreWithMapGeneral(MapScoreWorkspace *workspace,
                                      bool plot)
{
	CrystalPtr crystal = workspace->crystal;
	std::vector<AtomPtr> selected = workspace->selectAtoms;

	bool first = (workspace->segment == FFTPtr());

	if (first)
	{
		workspace->segment = prepareMapSegment(crystal, selected,
		                                       &workspace->basis,
		                                       &workspace->ave);
		
		workspace->constant = FFTPtr(new FFT(*workspace->segment));
		workspace->fcSegment = FFTPtr(new FFT(*workspace->segment));
		workspace->segment->createFFTWplan(1);
	}
	else
	{
		workspace->segment->copyFrom(workspace->constant);
	}
	
	/* Now we can add neighbouring atoms from the same Crystal
	* with impunity, and they won't go over the borderline already
	* established.*/
	if (first)
	{
		/* Make half-box measures */
		double xAng, yAng, zAng;
		xAng = workspace->segment->nx * workspace->segment->scales[0] / 2;
		yAng = workspace->segment->ny * workspace->segment->scales[1] / 2;
		zAng = workspace->segment->nz * workspace->segment->scales[2] / 2;
		
		if (!(workspace->flag & MapScoreFlagNoSurround))
		{
			workspace->extra = crystal->getAtomsInBox(workspace->ave, 
			                                          xAng, yAng, zAng);
			int count = 0;
			
			std::vector<AtomPtr> added;

			/* We want to add anything which is static to the 
			 * constant fraction. */
			for (size_t i = 0; i < workspace->extra.size(); i++)
			{
				AtomPtr anAtom = workspace->extra[i];

				if (std::find(selected.begin(), selected.end(), anAtom) 
				    != selected.end())
				{
					continue;
				}

				if (std::find(added.begin(), added.end(), anAtom) 
				    != added.end())
				{
					continue;
				}

				workspace->extra[i]->addToMap(workspace->constant, 
				                              workspace->basis,
				                              workspace->ave, 
				                              false, false, true);

				added.push_back(workspace->extra[i]);
				count++;
			}
			
			/* Copy this constant fraction into the segment */
			workspace->segment->copyFrom(workspace->constant);
		}
	}
	

	for (size_t i = 0; i < selected.size(); i++)
	{
		selected[i]->addToMap(workspace->segment, workspace->basis, 
		                      workspace->ave, false, false, true);
	}
	
	if (false && workspace->flag & MapScoreFlagReplaceWithObs)
	{
		FFTPtr fcSegment = workspace->segment;
		FFTPtr obsMap = crystal->getFFT();

		FFT::operation(obsMap, fcSegment, workspace->ave,
		               MapScoreTypeCopyToSmaller, NULL);
		
		return 0;
	}
	
	double score = 0;

	if (workspace->flag & MapScoreFlagSkipScore)
	{
		return 0;
	}
	
	/* In the middle of making calculation really quick?? */
	/* What the fuck did that comment mean? ^ */
	bool difference = (workspace->flag & MapScoreFlagDifference);
	ScoreType type = workspace->scoreType;
	
	if (workspace->flag & MapScoreFlagReplaceWithObs)
	{
		type = ScoreTypeAddDensity;
	}
	
	score = scoreFinalMap(crystal, workspace->segment, plot,
	                      type, workspace->ave, difference);

	return score;
}

double AtomGroup::scoreFinalValues(std::vector<double> xs,
                                   std::vector<double> ys,
                                   ScoreType scoreType,
                                   bool ignoreCutoff)
{
	double cutoff = MAP_VALUE_CUTOFF;
	
	if (ignoreCutoff) cutoff = -FLT_MAX;
	
	if (scoreType == ScoreTypeCorrel)
	{
		double correl = correlation(xs, ys, cutoff);
		return -correl;
	}
	else if (scoreType == ScoreTypeHappiness)
	{
		double happiness = happiness_coefficient(xs, ys);
		return -happiness;
	}
	else if (scoreType == ScoreTypeRFactor)
	{
		double rFactor = scaled_r_factor(xs, ys, cutoff);
		return rFactor;
	}
	else if (scoreType == ScoreTypeMultiply)
	{
		double mult = weightedMapScore(xs, ys);
		return -mult;
	}
	else if (scoreType == ScoreTypeAddDensity)
	{
		double sum = add_if_y_gt_zero(xs, ys);
		return sum;
	}
	else if (scoreType == ScoreTypeAddVoxels)
	{
		double sum = add_if_gt_zero(ys);
		return sum;
	}
	
	
	
	return 0;
}

void AtomGroup::plotCoordVals(std::vector<CoordVal> &vals, 
                              bool difference, double cutoff,
                              std::string filename)
{
	CSVPtr csv = CSVPtr(new CSV(6, "x", "y", "z", "fo", "fc", "mask"));

	for (size_t i = 0; i < vals.size(); i++)
	{
		double fo = vals[i].fo;
		double fc = vals[i].fc;
		double mask = 0;
		vec3 pos = make_vec3(0, 0, 0);

		#ifdef COORDVAL_FULL
		mask = vals[i].mask;
		pos = vals[i].pos;
		#endif

		if (!difference && fc < cutoff) continue;

		csv->addEntry(6, pos.x, pos.y, pos.z, fo, fc, mask);
	}

	csv->writeToFile(filename + ".csv");

	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = filename;
	plotMap["xHeader0"] = "fc";
	plotMap["yHeader0"] = "fo";
	plotMap["colour0"] = "black";

	plotMap["xTitle0"] = "Calc density";
	plotMap["yTitle0"] = "Obs density";
	plotMap["style0"] = "scatter";
	csv->plotPNG(plotMap);
}

double AtomGroup::scoreFinalMap(CrystalPtr crystal, FFTPtr segment,
                                bool plot, ScoreType scoreType,
                                vec3 ave, bool difference)
{
	double cutoff = MAP_VALUE_CUTOFF;
	mat3x3 real2Frac = crystal->getReal2Frac();

	/* Convert real2Frac to crystal coords to get correct segment
	* of the big real space map. */
	mat3x3_mult_vec(real2Frac, &ave);

	std::vector<double> xs, ys;
	std::vector<CoordVal> vals;

	FFTPtr map = crystal->getFFT();
	
	if (difference)
	{
		map = crystal->getDiFFT();
	}
	
	MapScoreType mapType = MapScoreTypeCorrel;

	if (scoreType == ScoreTypeAddDensity)
	{
		mapType = MapScoreTypeCopyToSmaller;
	}

	FFT::operation(map, segment, ave, mapType, &vals, true);

	/* For correlation calculations */
	for (size_t i = 0; i < vals.size(); i++)
	{
		xs.push_back(vals[i].fo);
		ys.push_back(vals[i].fc);
	}

	if (scoreType == ScoreTypeRFactor)
	{
		double scale = scale_factor_cutoff(xs, ys, cutoff);

		cutoff /= scale;	
		for (size_t i = 0; i < ys.size(); i++)
		{
			ys[i] /= scale;
			vals[i].fc /= scale;
		}
	}

	/* Debugging ... writes cc_score.csv and cc_score.png, csv can be
	* looked at with gnuplot quite nicely.*/

	if (plot)
	{
		plotCoordVals(vals, difference, cutoff, "cc_score");
	}


	/* Clear out the massive vectors */
	vals.clear();
	std::vector<CoordVal>().swap(vals);

	double score = scoreFinalValues(xs, ys, scoreType, difference);

	return score;
}

void AtomGroup::addProperties()
{
	addBoolProperty("been_tied", &_beenTied);
	addIntProperty("times_refined", &_timesRefined);

	for (size_t i = 0; i < atomCount(); i++)
	{
		addReference("atom", atom(i));
	}

}

void AtomGroup::addObject(ParserPtr object, std::string category)
{
	if (category == "atom")
	{
		AtomPtr atom = ToAtomPtr(object);
		addAtom(atom);
	} 
}

void AtomGroup::linkReference(ParserPtr object, std::string category)
{
	if (category == "atom")
	{
		AtomPtr atom = ToAtomPtr(object);
		addAtom(atom);
	}
}

std::vector<AtomPtr> AtomGroup::getHydrogenBonders()
{
	std::vector<AtomPtr> returns;
	for (int i = 0; i < atomCount(); i++)
	{
		bool add = atom(i)->canBeHydrogenBonder();

		if (add)
		{
			returns.push_back(atom(i));
		}
	}
	
	return returns;
}

void AtomGroup::refreshBondAngles()
{
	for (int i = 0; i < atomCount(); i++)
	{
		atom(i)->refreshBondAngles();
	}
}

vec3 AtomGroup::centroid()
{
	vec3 sum = empty_vec3();
	
	for (int i = 0; i < atomCount(); i++)
	{
		vec3 abs = atom(i)->getAbsolutePosition();
		vec3_add_to_vec3(&sum, abs);
	}
	
	vec3_mult(&sum, 1 / (double)atomCount());
	
	return sum;
}
