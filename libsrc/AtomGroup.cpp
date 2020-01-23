//
//  AtomGroup.cpp
//  vagabond
//
//  Created by Helen Ginn on 25/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#define MAP_VALUE_CUTOFF (0.005)

#include "Timer.h"
#include "Twist.h"
#include "ExplicitModel.h"
#include "AtomGroup.h"
#include <climits>
#include <algorithm>
#include "Atom.h"
#include "Anchor.h"
#include "Bond.h"
#include <sstream>
#include "Crystal.h"
#include <iomanip>
#include "CSV.h"
#include "maths.h"
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
	AtomList list;

	for (size_t i = 0; i < atomCount(); i++)
	{
		if ((atom(i)->getAtomName() == atomType)
		 && (atom(i)->getResidueNum() == resNum))
		{
			list.push_back(atom(i));
		}
	}

	return list;
}

AtomList AtomGroup::findAtomByNum(std::string atomType, int atomNum)
{
	AtomList list;

	for (size_t i = 0; i < atomCount(); i++)
	{
		if ((atom(i)->getAtomName() == atomType)
		 && (atom(i)->getAtomNum() == atomNum))
		{
			list.push_back(atom(i));
		}
	}

	return list;
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
	
	std::string add = "+" + group->getName();
	addToName(add);
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

		double val = atom(i)->posDisplacement(false, false);

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

void AtomGroup::propagateChangeExceptAnchor()
{
	for (size_t i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getModel()->isAnchor())
		{
			continue;
		}

		atom(i)->getModel()->propagateChange(0);
	}
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
	}

	for (size_t i = 0; i < atomCount(); i++)
	{
		if (!atom(i)) continue;

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
//		std::vector<AtomPtr>::iterator it;
//		it = std::lower_bound(_atoms.begin(), _atoms.end(), 
//		                      atom, Atom::greater);
//		_atoms.insert(it, atom);
		_atoms.push_back(atom);
		crystal->updateLargestNum(atom);
	}
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
					continue;
				}
				
				BondPtr bond = ToBondPtr(topAtom->getModel());

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

double AtomGroup::scoreWithMap(ScoreType scoreType, CrystalPtr crystal, 
                               std::string plot, unsigned int flags)
{
	OptionsPtr options = Options::getRuntimeOptions();
	DiffractionPtr data = options->getActiveData();
	if (!data || !atomCount())
	{
		return 0;
	}
	
	MapScoreWorkspace workspace;
	workspace.scoreType = scoreType;
	workspace.crystal = crystal;
	workspace.selectAtoms = shared_from_this();
	workspace.segment = FFTPtr();
	workspace.ave = empty_vec3();
	workspace.basis = make_mat3x3();
	workspace.flag = flags;
	workspace.filename = plot;

	return scoreWithMapGeneral(&workspace, plot != "");
}

FFTPtr AtomGroup::prepareMapSegment(CrystalPtr crystal,
                                    AtomGroupPtr selected,
                                    mat3x3 *basis, vec3 *ave)
{
	double maxDistance = 0;
	FFTPtr map = crystal->getFFT();
	mat3x3 real2Frac = crystal->getReal2Frac();

	vec3 sum = make_vec3(0, 0, 0);

	/* Find centroid of atom set */
	for (size_t i = 0; i < selected->atomCount(); i++)
	{
		selected->atom(i)->getModel()->refreshPositions();
		vec3 offset = selected->atom(i)->getModel()->getAbsolutePosition();
		sum = vec3_add_vec3(sum, offset);	
	}

	vec3_mult(&sum, 1 / (double)selected->atomCount());

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

	for (size_t i = 0; i < selected->atomCount(); i++)
	{
		/* Refresh absolute position */
		selected->atom(i)->getModel()->refreshPositions();
		vec3 offset = selected->atom(i)->getModel()->getAbsolutePosition();

		vec3 diff = vec3_subtract_vec3(offset, sum);

		ns[0] = std::max(fabs(diff.x), ns[0]);
		ns[1] = std::max(fabs(diff.y), ns[1]);
		ns[2] = std::max(fabs(diff.z), ns[2]);
	}

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

	vec3 half_box = make_vec3(nl[0] / 2, nl[1] / 2, nl[2] / 2);
//	vec3_subtract_from_vec3(ave, half_box);
	
	mat3x3_scale(basis, nl[0], nl[1], nl[2]);
	*basis = mat3x3_inverse(*basis);

	return segment;
}

void expand_limits(vec3 *min, vec3 *max, vec3 abs)
{
	if (abs.x < min->x) min->x = abs.x;
	if (abs.x > max->x) max->x = abs.x;
	if (abs.y < min->y) min->y = abs.y;
	if (abs.y > max->y) max->y = abs.y;
	if (abs.z < min->z) min->z = abs.z;
	if (abs.z > max->z) max->z = abs.z;
}

void AtomGroup::xyzLimits(vec3 *min, vec3 *max)
{
	*min = make_vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	*max = make_vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);

	if (atomCount() == 0)
	{
		*min = empty_vec3();
		*max = empty_vec3();
	}

	CrystalPtr crystal = Options::getActiveCrystal();
	for (int i = 0; i < atomCount(); i++)
	{
		vec3 abs = atom(i)->getAbsolutePosition();
		expand_limits(min, max, abs);

		if (atom(i)->getModel()->hasExplicitPositions())
		{
			for (int j = 0; j < crystal->getSampleNum(); j++)
			{
				ExplicitModelPtr m = atom(i)->getExplicitModel();
				vec3 abs = m->getFinalPositions()[j].start;

				expand_limits(min, max, abs);
			}
		}
	}
}

void AtomGroup::addToCubicMap(FFTPtr scratchFull, vec3 offset,
                              EleCache *cache)
{
	size_t nElements = totalElements();

	FFTPtr scratch = FFTPtr(new FFT(*scratchFull));
	scratch->wipe();
	scratch->takePlansFrom(scratchFull);
	
	FFTPtr eleFFT;
	if (cache == NULL)
	{
		eleFFT = FFTPtr(new FFT(*scratchFull));
		eleFFT->wipe();
	}

	std::vector<AtomPtr> after;
	size_t totalElectrons = 0;

	mat3x3 new_basis = scratchFull->getReal2Frac();

	for (int i = 0; i < nElements; i++)
	{
		ElementPtr ele = element(i);
		
		if (cache != NULL && cache->count(ele))
		{
			eleFFT = (*cache)[ele];
		}
		else
		{
			if (!eleFFT)
			{
				eleFFT = FFTPtr(new FFT(*scratchFull));
			}

			eleFFT->wipe();
			ele->populateFFT(new_basis, eleFFT);
			
			if (cache != NULL)
			{
				(*cache)[ele] = eleFFT;
			}
		}

		double elementElectrons = 0;

		for (int j = 0; j < atomCount(); j++)
		{
			if (atom(j)->getElement() == ele)
			{
				atom(j)->addDirectlyToMap(scratch, new_basis, offset);
				double eCount = ele->electronCount() *
				atom(j)->getModel()->getEffectiveOccupancy();
				
				elementElectrons += eCount;
			}
		}
		
		totalElectrons += elementElectrons;
		
		scratch->fft(1);
		FFT::multiply(scratch, eleFFT);
		scratch->fft(-1);
		scratch->setTotal(elementElectrons);
		FFT::addSimple(scratchFull, scratch);
		scratch->wipe();
	}
	
	FFT::addSimple(scratchFull, scratch);
}

vec3 AtomGroup::prepareCubicMap(FFTPtr *scratchFull, vec3 *offset, 
                                vec3 min, vec3 max, double buff)
{
	double maxDStar = Options::getRuntimeOptions()->getActiveCrystalDStar();
	double cubeDim = 1.0 / (2.0 * maxDStar);
	
	const double maxD = 0.8;
	if (cubeDim > maxD)
	{
		cubeDim = maxD;
	}
	
	CrystalPtr crystal = Options::getActiveCrystal();
	

	/* 2 Angstroms buffer region on either side of the protein */
	vec3 buffer = make_vec3(buff, buff, buff);
	vec3_subtract_from_vec3(&min, buffer);

	/* Modify the offset to include buffer region */
	vec3_add_to_vec3(offset, min);

	vec3_add_to_vec3(&max, buffer);
	vec3 limits = vec3_subtract_vec3(max, min);

	/* find out how many voxels in each dimension to cover entire atom group
	 * with a regular grid */
	vec3 extent = limits;
	vec3_mult(&extent, 1 / cubeDim);
	extent.x = (int)lrint(extent.x);
	extent.y = (int)lrint(extent.y);
	extent.z = (int)lrint(extent.z);
	
	(*scratchFull) = FFTPtr(new FFT());
	(*scratchFull)->create(extent.x, extent.y, extent.z);
	(*scratchFull)->setScales(cubeDim);
	(*scratchFull)->createFFTWplan(1);
	
	return empty_vec3();
}

vec3 finalOffset(vec3 offset, FFTPtr cube,
                 mat3x3 real2frac, double buff = BUFFER_REGION)
{
	/* Add half box before addition to the FFT of arbitrary voxel size */
	double scale = cube->getScale(0);
	vec3 half_box = make_vec3(cube->nx, cube->ny, cube->nz);
	vec3_mult(&half_box, scale * 0.5);
	
	vec3_add_to_vec3(&offset, half_box);
	vec3 buffer = make_vec3(buff, buff, buff);
	vec3_subtract_from_vec3(&offset, buffer);
	mat3x3_mult_vec(real2frac, &offset);

	return offset;
}

void AtomGroup::addToMap(FFTPtr fft, mat3x3 real2frac, vec3 offset,
                         EleCache *cache)
{
	size_t nElements = totalElements();
	
	if (nElements == 0)
	{
		return;
	}
	
	FFTPtr scratchFull;

	double buffer = BUFFER_REGION;
	vec3 min, max;
	xyzLimits(&min, &max);
	prepareCubicMap(&scratchFull, &offset, min, max, buffer);
	addToCubicMap(scratchFull, offset);
	scratchFull->fft(1);
	double b = Options::getActiveCrystal()->getRealBFactor();
	Distributor::bFactorDistribute(scratchFull, b);
	scratchFull->fft(-1);

	vec3 addvec = finalOffset(min, scratchFull, real2frac, buffer);
	
	FFT::operation(fft, scratchFull, addvec, MapScoreTypeNone, 
	               NULL, false, true);

	
}

double AtomGroup::scoreWithMapGeneral(MapScoreWorkspace *workspace,
                                      bool plot)
{
	if (plot && workspace->filename == "")
	{
		workspace->filename = "cc_score.png";
	}

	if (!workspace->selectAtoms->atomCount())
	{
		return 0;
	}

	CrystalPtr crystal = workspace->crystal;
	AtomGroupPtr selected = workspace->selectAtoms;

	bool first = (workspace->segment == FFTPtr());

	vec3 min, max; 
	selected->xyzLimits(&min, &max);

	if (first)
	{
		workspace->ave = empty_vec3();
		workspace->eleCache.clear();
		selected->prepareCubicMap(&workspace->segment,
		                          &workspace->ave, min, max);
		workspace->basis = workspace->segment->getReal2Frac();

		workspace->constant = FFTPtr(new FFT(*workspace->segment));
		workspace->constant->takePlansFrom(workspace->segment);
	}
	else
	{
//		workspace->segment->copyFrom(workspace->constant);
			workspace->segment->wipe();
	}
	
	/* Now we can add neighbouring atoms from the same Crystal
	 * once we have established they do not go over the border. */
	if (first)
	{
		/* Make half-box measures */
		vec3 minmax = vec3_add_vec3(min, max);
		vec3_mult(&minmax, 0.5);
		vec3 half = vec3_subtract_vec3(max, min);
		vec3_mult(&half, 0.5);
		
		if (!(workspace->flag & MapScoreFlagNoSurround))
		{
			workspace->extra = crystal->getAtomsInBox(minmax, half.x, 
			                                          half.y, half.z);

			AtomGroupPtr added = AtomGroupPtr(new AtomGroup());

			/* We want to add anything which is static to the 
			 * constant fraction. */
			for (size_t i = 0; i < workspace->extra->atomCount(); i++)
			{
				AtomPtr anAtom = workspace->extra->atom(i);

				if (workspace->selectAtoms->hasAtom(anAtom) ||
				    added->hasAtom(anAtom))
				{
					continue;
				}

				added->addAtom(anAtom);
			}
			
			/* Add the acceptable atoms to the map */
			added->addToCubicMap(workspace->constant, 
			                     workspace->ave,
			                     &workspace->eleCache);
			
			/* Copy this constant fraction into the segment */
//			workspace->segment->copyFrom(workspace->constant);
		}
	}
	
	workspace->selectAtoms->addToCubicMap(workspace->segment,
	                                      workspace->ave,
	                                      &workspace->eleCache);
	
	double score = 0;

	if (workspace->flag & MapScoreFlagSkipScore)
	{
		return 0;
	}
	
	/* In the middle of making calculation really quick?? */
	/* What the fuck did that comment mean? ^ */
	bool difference = (workspace->flag & MapScoreFlagDifference);
	ScoreType type = workspace->scoreType;
	
	/* Convert average for comparing */
	if (first)
	{
		mat3x3 r2f = crystal->getReal2Frac();
		workspace->working_ave = finalOffset(min, workspace->segment, r2f);
	}

	score = scoreFinalMap(workspace, plot);

	return score;
}

double AtomGroup::scoreFinalMap(MapScoreWorkspace *ws, bool plot)
{
	/* Convert real2Frac to crystal coords to get correct segment
	* of the big real space map. */
	CrystalPtr crystal = ws->crystal;
	vec3 ave = ws->working_ave;
	
	double cutoff = MAP_VALUE_CUTOFF;

	std::vector<double> xs, ys, weights;
	std::vector<CoordVal> vals;

	FFTPtr map = ws->crystal->getFFT();
	bool difference = (ws->flag & MapScoreFlagDifference);
	
	if (difference)
	{
		map = crystal->getDiFFT();
	}
	
	MapScoreType mapType = MapScoreTypeCorrel;
	
	if (ws->scoreType == ScoreTypeCopyToSmaller)
	{
		mapType = MapScoreTypeCopyToSmaller;
	}
	
	/* Actual variables to compare, correlation calculations */
	ws->segment->copyRealToImaginary();
	
	ws->vals.clear();
	FFT::addSimple(ws->segment, ws->constant); /* add to real */
	FFT::operation(map, ws->segment, ave, mapType, &ws->vals);

	for (size_t i = 0; i < ws->vals.size() && 
	     ws->scoreType != ScoreTypeCorrel; i++)
	{
		if (ws->vals[i].weight > 1e-6)
		{
			ys.push_back(ws->vals[i].fc);
			xs.push_back(ws->vals[i].fo);
			double weight = (ws->vals[i].weight);
			weights.push_back(weight);
		}
	}

	/* Debugging ... writes cc_score.csv and cc_score.png, csv can be
	* looked at with gnuplot quite nicely.*/

	if (plot)
	{
		plotCoordVals(ws->vals, difference, cutoff, ws->filename);
	}
	
	if (ws->scoreType == ScoreTypeCorrel)
	{
		double correl = correlation(ws->vals);
		return -correl;
	}

	if (ws->scoreType == ScoreTypeRFactor)
	{
		double rfactor = weighted_r_factor(ws->vals);
		return rfactor;
	}

	double score = scoreFinalValues(xs, ys, weights, ws->scoreType, ws->flag);

	return score;
}


double AtomGroup::scoreFinalValues(std::vector<double> &xs,
                                   std::vector<double> &ys,
                                   std::vector<double> &weights,
                                   ScoreType scoreType,
                                   unsigned int flags)
{
	bool difference = (flags & MapScoreFlagDifference);
	double cutoff = MAP_VALUE_CUTOFF;
	
//	if (difference) cutoff = -FLT_MAX;
	cutoff = -FLT_MAX;
	
	if (scoreType == ScoreTypeCorrel)
	{
		double correl = correlation(xs, ys, cutoff, &weights);
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
	if (filename.length() == 0)
	{
		filename = "cc_score";
	}

	CSVPtr csv = CSVPtr(new CSV(5, "x", "y", "z", "fo", "fc"));

	for (size_t i = 0; i < vals.size(); i++)
	{
		double fo = vals[i].fo;
		double fc = vals[i].fc;
		vec3 pos = make_vec3(0, 0, 0);

		#ifdef COORDVAL_FULL
		pos = vals[i].pos;
		#endif

		if (!difference && fc < cutoff) continue;

		csv->addEntry(5, pos.x, pos.y, pos.z, fo, fc);
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

void AtomGroup::addProperties()
{
	addBoolProperty("been_tied", &_beenTied);
	addIntProperty("times_refined", &_timesRefined);
	addStringProperty("name", &_name);

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

void AtomGroup::linkReference(BaseParserPtr object, std::string category)
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

void AtomGroup::makeBackboneTwists(ExplicitModelPtr applied)
{
	for (int i = 0; i < atomCount(); i++)
	{
		AtomPtr a = atom(i);
		
		if (!a->isBackbone() && !a->isBackboneAndSidechain())
		{
			continue;
		}
		
		if (!a->getModel()->isBond())
		{
			continue;
		}
		
		BondPtr b = ToBondPtr(a->getModel());
		
		if (b->isFixed())
		{
			continue;
		}

		TwistPtr twist = TwistPtr(new Twist());
		twist->setBond(b);
		twist->addToAppliedModel(applied);
	}
}

void AtomGroup::boundingMonomers(int *begin, int *end)
{
	*begin = INT_MAX;
	*end = -INT_MAX;

	for (int i = 0; i < atomCount(); i++)
	{
		AtomPtr a = atom(i);
		int resi = a->getResidueNum();
		
		if (*begin > resi)
		{
			*begin = resi;
		}
		
		if (*end < resi)
		{
			*end = resi;
		}
	}
}

size_t AtomGroup::totalElements()
{
	_elements.clear();

	for (int i = 0; i < atomCount(); i++)
	{
		ElementPtr ele = atom(i)->getElement();
		if (std::find(_elements.begin(), _elements.end(), ele) 
		    == _elements.end())
		{
			_elements.push_back(ele);
		}
	}

	return _elements.size();
}

double AtomGroup::recalculatePositions(void *obj)
{
	static_cast<AtomGroup *>(obj)->refreshPositions();

	return 0;
}

