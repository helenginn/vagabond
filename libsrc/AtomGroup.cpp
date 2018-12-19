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
#include "Anchor.h"
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

AtomGroupPtr AtomGroup::subGroupForConf(int conf)
{
	AtomGroupPtr group = AtomGroupPtr(new AtomGroup());
	std::string confID = conformer(conf);
	
	for (size_t i = 0; i < atomCount(); i++)
	{
		if (_atoms[i]->getAlternativeConformer() == confID)
		{
			group->addAtom(_atoms[i]);
		}
	}

	return group;

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

	return conformerMap().size();
}

int AtomGroup::conformer(std::string conf)
{
	std::map<std::string, size_t> conformerList = conformerMap();
	
	if (!conformerList.count(conf))
	{
		return 0;
	}
	else
	{
		return conformerList[conf];
	}
}

std::string AtomGroup::conformer(size_t i)
{
	if (i > conformerCount()) return "";

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

		stream << atom(i)->anisouPDBLine(crystal);
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
		AtomPtr atom = list[i];
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

AtomList AtomGroup::beyondGroupAtoms(bool just_bottom)
{
	AtomList list = topLevelAtoms();
	AtomList bottom;
	
	for (int i = 0; i < list.size(); i++)
	{
		AtomPtr a = list[i];
		
		if (!a->getModel()->isBond())
		{
			continue;
		}
		
		AtomPtr last = a;

		while (hasAtom(a))
		{
			BondPtr b = ToBondPtr(a->getModel());

			if (!(b->downstreamBondGroupCount() && b->downstreamBondCount(0)))
			{
				a = AtomPtr();
				break;
			}
			
			b = ToBondPtr(b->downstreamBond(0, 0));
			last = a;
			a = b->getMinor();
		}
		
		if (a && !just_bottom)
		{
			bottom.push_back(a);
		}
		else if (last && just_bottom)
		{
			bottom.push_back(last);
		}
	}
	
	return bottom;
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

		while (true)
		{
			if (topAtom->getModel()->isBond() &&
			    topAtom->getAlternativeConformer() == conf)
			{
				break;
			}
			
			j++;

			if (j >= atomCount())
			{
				goto giveup;
			}

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
			
			if (bond->getMajor()->getModel()->isAnchor())
			{
				list.push_back(topAtom);
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
	if (!anAtom) return false;

	for (size_t i = 0; i < atomCount(); i++)
	{
		if (atom(i) == anAtom)
		{
			return true;
		}
	}

	return false;
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

void AtomGroup::saveAtomPositions()
{
	for (int i = 0; i < atomCount(); i++)
	{
		atom(i)->saveInitialPosition();
	}
}

void AtomGroup::removeAtom(AtomPtr atom)
{
	if (!atom)
	{
		return;
	}

	std::vector<AtomPtr>::iterator it;
	it = std::find(_atoms.begin(), _atoms.end(), atom);

	if (it != _atoms.end())
	{
		_atoms.erase(it);
	}
}

void AtomGroup::addAtom(AtomPtr atom)
{
	if (!atom)
	{
		return;
	}

	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();

	std::vector<AtomPtr>::iterator it;
	it = std::find(_atoms.begin(), _atoms.end(), atom);

	if (it == _atoms.end())
	{
		_atoms.push_back(atom);
		crystal->updateLargestNum(atom);
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
	if (!atomCount())
	{
		return;
	}
	
	AtomList topAtoms = topLevelAtoms();
	
	_timesRefined++;

	ScoreType scoreType = ScoreTypeModelPos;
	int maxTries = 0;

	switch (rType) {
		case RefinementModelPos:
		scoreType = ScoreTypeModelPos;
		maxTries = 60;
		break;

		case RefinementSavedPos:
		scoreType = ScoreTypeSavedPos;
		maxTries = 200;
		break;

		case RefinementCentroid:
		scoreType = ScoreTypeCentroid;
		maxTries = 60;
		break;

		case RefinementCrude:
		scoreType = ScoreTypeCorrel;
		maxTries = 1;
		break;

		case RefinementFine:
		scoreType = ScoreTypeCorrel;
		maxTries = 2;
		break;

		case RefinementMouse:
		scoreType = ScoreTypeMouse;
		maxTries = 1;
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
	
	if (hasParameter(ParamOptionMaxTries))
	{
		maxTries = getParameter(ParamOptionMaxTries);
	}
	
	while (topAtoms.size() > 0)
	{
		BondPtr topBond;
		int count = 0;

		bool changed = true;

		while (changed && count < maxTries)
		{
			setupNelderMead();
			setCrystal(target);
			setCycles(16);
			setScoreType(scoreType);

			for (int i = 0; i < topAtoms.size(); i++)
			{
				AtomPtr topAtom = topAtoms[i];

				if (!topAtom->getModel()->isBond())
				{
					if (rType == RefinementSavedPos &&
					    topAtom->getModel()->isAnchor() && false)
					{
						AnchorPtr anch = ToAnchorPtr(topAtom->getModel());
						addAnchorPosition(anch, 0.02, 0.002);
						
						if (hasParameter(ParamOptionTTN))
						{
							topAtom = anch->getNAtom();
						}
						else if (hasParameter(ParamOptionTTC))
						{
							topAtom = anch->getCAtom();
						}
						else
						{
							continue;
						}
					}
					else
					{
						continue;
					}
				}
				
				BondPtr bond = ToBondPtr(topAtom->getModel());
				backwards = isBackwards(bond);

				if (i == 0)
				{
					setJobName("torsion_" +  bond->shortDesc());
				}

				if (rType != RefinementFine && 
				    rType != RefinementMouse &&
				    rType != RefinementCrude)
				{
					topBond = setupThoroughSet(bond, false);
				}
				else
				{
					topBond = setupThoroughSet(bond);
				}
			}

			for (size_t l = 0; l < _includeForRefine.size(); l++)
			{
				addSampledAtoms(_includeForRefine[l]);
			}

			if (rType == RefinementModelPos 
			    || rType == RefinementSavedPos 
			    || rType == RefinementFine 
			    || rType == RefinementCrude 
			    || rType == RefinementModelRMSDZero
			    || rType == RefinementRMSDZero
			    || rType == RefinementMouse)
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
		
		if (backwards)
		{
//			std::cout << "I am backwards!" << std::endl;
			
			ModelPtr model = topBond->getParentModel();
			if (!model->isBond() && !model->isAnchor())
			{
				break;
			}

			topAtom = ToBondPtr(model)->getMajor();
		}
		
		if (rType == RefinementSavedPos && !backwards)
		{
//			std::cout << "I am not backwards" << std::endl;
		}
		
		if (!topAtom)
		{
			break;
		}
		
		topAtoms = findAtoms(topAtom->getAtomName(), topAtom->getResidueNum());
		
		if (rType == RefinementCrude)
		{
			return;
		}
	}

	_includeForRefine.clear();
	refreshPositions(true);
}

double AtomGroup::scoreWithMap(ScoreType scoreType, CrystalPtr crystal, bool plot, unsigned int flags)
{
	OptionsPtr options = Options::getRuntimeOptions();
	DiffractionPtr data = options->getActiveData();
	if (!data || !atomCount())
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
	if (!workspace->selectAtoms.size())
	{
		return 0;
	}
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
	
	score = scoreFinalMap(crystal, workspace->segment, plot,
	                      type, workspace->ave, workspace->flag);

	return score;
}

double AtomGroup::scoreFinalValues(std::vector<double> xs,
                                   std::vector<double> ys,
                                   ScoreType scoreType,
                                   unsigned int flags)
{
	bool difference = (flags & MapScoreFlagDifference);
	double cutoff = MAP_VALUE_CUTOFF;
	
	if (difference) cutoff = -FLT_MAX;
	
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
		int option = 0;
		
		if (flags & MapScoreFlagPosOnly)
		{
			option = 1;
		}
		else if (flags & MapScoreFlagNegOnly)
		{
			option = -1;
		}
		
		double sum = add_x_if_y(xs, ys, option);
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
                                vec3 ave, unsigned int flags)
{
	
	double cutoff = MAP_VALUE_CUTOFF;
	mat3x3 real2Frac = crystal->getReal2Frac();

	/* Convert real2Frac to crystal coords to get correct segment
	* of the big real space map. */
	mat3x3_mult_vec(real2Frac, &ave);

	std::vector<double> xs, ys;
	std::vector<CoordVal> vals;

	FFTPtr map = crystal->getFFT();
	bool difference = (flags & MapScoreFlagDifference);
	
	if (difference)
	{
		map = crystal->getDiFFT();
	}
	
	MapScoreType mapType = MapScoreTypeCorrel;
	
	if (scoreType == ScoreTypeCopyToSmaller)
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

	double score = scoreFinalValues(xs, ys, scoreType, flags);

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
