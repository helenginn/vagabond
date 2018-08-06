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
#include "Element.h"
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

AtomPtr AtomGroup::findAtom(std::string atomType, std::string confID)
{
	AtomList atoms = findAtoms(atomType);

	for (size_t i = 0; i < atoms.size(); i++)
	{
		if (atoms[i].expired())
		{
			continue;
		}

		if (atoms[i].lock()->getAlternativeConformer() == confID)
		{
			return atoms[i].lock();
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
		total += atom(i)->getElement()->electronCount();
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
		/* Get the total number of conformers to worry about */
		std::vector<BondSample> *samples = atom(0)->getModel()->getManyPositions();

		numConf = samples->size();

		for (size_t i = 0; i < atomCount(); i++)
		{
			if (atom(i)->getWeighting() <= 0)
			{
				continue;
			}

			if (atom(i)->getModel()->hasExplicitPositions())
			{
				stream << atom(i)->getPDBContribution(conformer);
			}
			else
			{
				stream << atom(i)->averagePDBContribution(false, false);
			}
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
		if (atom(i)->getElement()->electronCount() <= 1)
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
		if (atom(i)->getElement()->electronCount() <= 1)
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
		atom(i)->getModel()->getFinalPositions();
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
		double e = atom(i)->getElement()->electronCount();
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

	for (size_t i = 0; i < 1; i++)
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

			topAtom = bond->getMajor();
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
		vec3 tmp = atom(i)->getAsymUnitPosition(crystal);
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
	int bondNum = 4;
	double degrees = 0;

	switch (rType) {
		case RefinementModelPos:
		scoreType = ScoreTypeModelPos;
		maxTries = 60;
		degrees = 4;
		break;

		case RefinementFine:
		scoreType = ScoreTypeCorrel;
		maxTries = 3;
		degrees = 4;
		bondNum = 4;
		refineAngles = false;
		break;

		case RefinementModelRMSDZero:
		scoreType = ScoreTypeModelRMSDZero;
		maxTries = 10;
		degrees = 4;
		break;

		case RefinementRMSDZero:
		scoreType = ScoreTypeRMSDZero;
		maxTries = 10;
		degrees = 4;
		break;

		default:
		shout_at_helen("Unimplemented refinement option?");
		break;
	}

	if (refineAngles)
	{
		bondNum = 3;
	}
	
	if (hasParameter(ParamOptionTorsion) && paramCount() < 2 &&
	    rType == RefinementFine)
	{
		maxTries = 20;
	}

	for (size_t n = 0; n < topAtoms.size(); n++)
	{
		AtomPtr topAtom = topAtoms[n].lock();

//		if (n > 0) std::cout << "'" << std::flush;

		while (hasAtom(topAtom))
		{
			if (!topAtom->getModel()->isBond())
			{
				break;
			}

			BondPtr bond = ToBondPtr(topAtom->getModel());

			int groups = bond->downstreamAtomGroupCount();

			if (!groups)
			{
				break;
			}

			if (!bond->isRefinable())
			{
				break;
			}

			int count = 0;

			BondPtr topBond;

			for (int k = 0; k < 1; /*bond->downstreamAtomGroupCount();*/ k++)
			{
				bool changed = true;
				bool addFlex = (rType == RefinementFine);

				while (changed && count < maxTries)
				{
					bond->setActiveGroup(k);
					setupNelderMead();
					setCrystal(target);
					setCycles(16);

					if (rType != RefinementFine)
					{
						topBond = setupTorsionSet(bond, k, bondNum,
						                          deg2rad(degrees), 
						                          deg2rad(0.04),
						                          refineAngles, addFlex);
					}
					else
					{
						topBond = setupThoroughSet(bond, bondNum,
						                           deg2rad(degrees), 
						                           deg2rad(0.04),
						                           refineAngles, addFlex);

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

					setJobName("torsion_" +  bond->shortDesc());
					changed = sample();
					count++;
				}

				if (!topBond)
				{
					topAtom = AtomPtr();
					continue;
				}

				topAtom = topBond->getMinor();
			}
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
		selected[i]->getModel()->getFinalPositions();
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
		selected[i]->getModel()->getFinalPositions();
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
		nl[i] = (2 * (ns[i] + buffer));
	}
	
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

	/* but the basis of the map workspace should be a mini-unit cell */
	double toReal[3];

	for (int i = 0; i < 3; i++)
	{
		toReal[i] = 1 / (segment->scales[i] * nl[i]);
	}

	mat3x3_scale(basis, nl[0], nl[1], nl[2]);
	*basis = mat3x3_inverse(*basis);
	

	return segment;
}

double AtomGroup::addAtomsQuickly(FFTPtr segment, std::vector<AtomPtr> selected, 
                                  mat3x3 basis, vec3 ave)
{
	std::vector<ElementPtr> elements = Element::elementList(selected);

	std::vector<AtomPtr> traditional;
	int allElementElectrons = 0;

	Timer tReal("real space");
	Timer tFFT("fft");

	FFTPtr tmpSegment = FFTPtr(new FFT(*segment));
	tmpSegment->setAll(0);
	tmpSegment->setupMask();
	tmpSegment->avoidWriteToMaskZero(false);

	/* For each element, make a new map and add all real space atom
	*  positions in one fell swoop, if possible */
	for (size_t i = 0; i < elements.size(); i++)
	{
		int totalElectrons = 0;
		FFTPtr elesegment = FFTPtr(new FFT(*segment));

		tReal.start();
		for (size_t j = 0; j < selected.size(); j++)
		{
			if (selected[j]->getElement() != elements[i])
			{
				continue;
			}

			ModelPtr model = selected[j]->getModel();
			model->getFinalPositions();

			/** This needs to be added in the old way at the end. */
			if (!model->hasExplicitPositions())
			{
				traditional.push_back(selected[j]);
				continue;
			}

			totalElectrons += elements[i]->electronCount();
			vec3 pos = selected[j]->getAbsolutePosition();
			pos = vec3_subtract_vec3(pos, ave);

			model->addRealSpacePositions(elesegment, pos);
		}
		tReal.stop();

		/* Must include models which are not explicit too */
		allElementElectrons += elements[i]->electronCount() * selected.size();

		elesegment->createFFTWplan(1);

		tFFT.start();
		elesegment->fft(1);
		tFFT.stop();

		FFTPtr elementDist = elements[i]->getDistribution(false,
		                                                  elesegment->nx);
		FFT::multiply(elesegment, elementDist);
		tFFT.start();
		elesegment->fft(-1);
		tFFT.stop();

		/* But this total electron count for a given element should not */
		elesegment->setTotal(totalElectrons * 10e4);

		FFT::addSimple(tmpSegment, elesegment);
	}

	//	tReal.report();
	//	tFFT.report();

	for (size_t i = 0; i < traditional.size(); i++)
	{
		traditional[i]->addToMap(tmpSegment, basis, ave);
	}

	tmpSegment->setTotal(allElementElectrons * 10e2);

	FFT::addSimple(segment, tmpSegment);
}

double AtomGroup::scoreWithMapQuick(ScoreType scoreType, CrystalPtr crystal,
                                    bool plot, std::vector<AtomPtr> selected)
{
	mat3x3 basis;
	vec3 ave;
	FFTPtr segment = prepareMapSegment(crystal, selected, &basis, &ave);

	addAtomsQuickly(segment, selected, basis, ave);

	double cutoff = MAP_VALUE_CUTOFF;
	segment->aboveValueToMask(cutoff);
	segment->avoidWriteToMaskZero();

	std::vector<AtomPtr> extra; 

	/* Neighbours */
	extra = crystal->getCloseAtoms(selected, 1.5);
	addAtomsQuickly(segment, extra, basis, ave);

	double score = scoreFinalMap(crystal, segment, plot, scoreType, ave);
	return score;
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

			workspace->segment->copyFrom(workspace->constant);
		}
	}

	for (size_t i = 0; i < selected.size(); i++)
	{
		selected[i]->addToMap(workspace->segment, workspace->basis, 
		                      workspace->ave, false, false, true);
	}
	
	if (workspace->flag & MapScoreFlagSubtractFc)
	{
		if (first)
		{
			FFTPtr fcSegment = workspace->fcSegment;
			FFTPtr map = crystal->getCalculatedMap();
			FFT::score(map, fcSegment, workspace->ave,
			           NULL, MapScoreTypeCopyToSmaller);
			fcSegment->scaleToFFT(workspace->segment);
			
			std::vector<CoordVal> vals;
			FFT::score(workspace->segment, fcSegment, empty_vec3(),
			           &vals, MapScoreTypeCorrel);
			plotCoordVals(vals, 0, 0, "cc_fcs");

			
			/* Get ready for subtraction! */
			fcSegment->multiplyAll(-1);
		}

		/* Remove the last calculated map from the segment, which
		 * is already negative. */
		FFT::addSimple(workspace->segment, workspace->fcSegment);	
	}
	
	double score = 0;

	/* In the middle of making calculation really quick?? */
	bool difference = (workspace->flag & MapScoreFlagDifference);
	
	score = scoreFinalMap(crystal, workspace->segment, plot,
	                      workspace->scoreType, workspace->ave,
	                      difference);

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

	plotMap["xTitle0"] = "Fc density";
	plotMap["yTitle0"] = "Fo density";
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

	FFT::operation(map, segment, ave, MapScoreTypeCorrel, &vals, true);

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
