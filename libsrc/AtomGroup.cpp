//
//  AtomGroup.cpp
//  vagabond
//
//  Created by Helen Ginn on 25/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#define MAP_VALUE_CUTOFF (0.005)

#include <hcsrc/Timer.h>
#include "vagmaths.h"
#include "Twist.h"
#include "Hydrogenator.h"
#include "Monomer.h"
#include "Element.h"
#include "ExplicitModel.h"
#include "AtomGroup.h"
#include <climits>
#include <algorithm>
#include "Blob.h"
#include "Atom.h"
#include "SymAtom.h"
#include "Anchor.h"
#include "Sponge.h"
#include "Bond.h"
#include <sstream>
#include "Crystal.h"
#include <iomanip>
#include "CSV.h"
#include "Bucket.h"
#include "Diffraction.h"
#include "Shouter.h"
#include "libccp4/ccp4_spg.h"
#include "Options.h"
#include "FFT.h"
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

AtomList AtomGroup::findAtoms(std::string atomType, int resNum,
                              std::string chainID)
{
	AtomList list;

	for (size_t i = 0; i < atomCount(); i++)
	{
		if ((atom(i)->getAtomName() == atomType
		 && (atom(i)->getResidueNum() == resNum) || resNum == INT_MAX))
		{
			if (chainID.length())
			{
				std::string chain = atom(i)->getMolecule()->getChainID();
				
				if (chain.substr(0, chainID.length()) != chainID)
				{
					continue;
				}

			}

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

void AtomGroup::addAtomsFrom(std::vector<AtomPtr> group)
{
	for (size_t i = 0; i < group.size(); i++)
	{
		addAtom(group[i]);	
	}
}

void AtomGroup::addAtomsFrom(AtomGroupPtr group)
{
	addAtomsFrom(&*group);
}

void AtomGroup::addAtomsFrom(AtomGroup *group)
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
		for (size_t i = 0; i < atomCount(); i++)
		{
			if (!atom(i)->getModel()->hasExplicitPositions())
			{
				stream << atom(i)->averagePDBContribution(false, false);
				continue;
			}

			if (atom(i)->getWeighting() <= 0)
			{
				continue;
			}

			if (conformer >= 0)
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
	_padding = 0;
	_t1 = new Timer("adding atoms");
	_t2 = new Timer("internal FFT addition");
	_t3 = new Timer("FFT transformation");
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
		atom(i)->getModel()->propagateChange(0);
	}

	for (size_t i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getMonomer() && 
		    atom(i)->getMonomer()->getIdentifier() == "pro")
		{
			if (atom(i)->getElement()->getSymbol() == "H")
			{
				if (atom(i)->getModel()->isBond())
				{
					BondPtr b = ToBondPtr(atom(i)->getModel());
					Hydrogenator::adjustProlineHydrogens(b, true);
				}
			}
		}

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

		/* Find an appropriate starting atom */
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

std::vector<AtomPtr> AtomGroup::getCloseAtoms(std::vector<AtomPtr> atoms, 
                                              double tol)
{
	std::vector<AtomPtr> extra;
	
	for (size_t i = 0; i < atoms.size(); i++)
	{
		std::vector<AtomPtr> clAtoms = getCloseAtoms(atoms[i], tol);
		
		for (int j = 0; j < clAtoms.size(); j++)
		{
			if (std::find(atoms.begin(), atoms.end(), clAtoms[j]) !=
			    atoms.end())
			{
				continue;
			}

			if (std::find(extra.begin(), extra.end(), clAtoms[j]) !=
			    extra.end())
			{
				continue;
			}

			extra.push_back(clAtoms[j]);
		}
	}
	
	return extra;
}

AtomGroupPtr AtomGroup::getAtomsInBox(vec3 target, double tolx,
                                      double toly, double tolz,
                                      bool addSyms)
{
	AtomGroupPtr atoms = AtomGroupPtr(new AtomGroup());
	
	int symopTry = 1;
	CrystalPtr crystal = DynCrystalPtr(getTopParser());
	
	if (!crystal) addSyms = false;

	mat3x3 real2Frac = crystal->getReal2Frac();
	mat3x3 frac2Real = crystal->getFrac2Real();
	vec3 frac = mat3x3_mult_vec(real2Frac, target);

	if (addSyms)
	{
		symopTry = crystal->symOpCount();
	}

	for (int j = 0; j < atomCount(); j++)
	{
		for (int i = 0; i < symopTry; i++)
		{
			AtomPtr anAtom = atom(j);
			vec3 pos = anAtom->getAbsolutePosition();
			vec3 next = anAtom->getSymRelatedPosition(i, pos);
			mat3x3_mult_vec(real2Frac, &next);
			vec3 copy_next = next;
			
			while (frac.x - next.x > 1) { next.x += 1; }
			while (frac.x - next.x <= -1) { next.x -= 1; }

			while (frac.y - next.y > 1)   { next.y += 1; }
			while (frac.y - next.y <= -1) { next.y -= 1; }

			while (frac.z - next.z > 1)   { next.z += 1; }
			while (frac.z - next.z <= -1) { next.z -= 1; }
			
			vec3 shift = vec3_subtract_vec3(next, copy_next);
			pos = mat3x3_mult_vec(frac2Real, next);

			vec3 diff = vec3_subtract_vec3(pos, target);

			if (fabs(diff.x) > tolx || fabs(diff.y) > toly
			    || fabs(diff.z) > tolz)
			{
				continue;
			}

			if (i == 0)
			{
				atoms->addAtom(anAtom);
			}
			else
			{
				SymAtomPtr satom = SymAtomPtr(new SymAtom(*anAtom));
				satom->setUnitCellShift(shift);
				satom->setSymop(i);
				atoms->addAtom(satom);
			}
		}
	}

	return atoms;
}

std::vector<AtomPtr> AtomGroup::getCloseAtoms(AtomPtr one, double tol, 
                                              bool cache)
{
	std::vector<AtomPtr> atoms;
	std::vector<AtomPtr> *searchAtoms = &_atoms;
	std::vector<AtomPtr> *atomListPtr = &atoms;
	
	if (cache)
	{
		atomListPtr = &_closeishAtoms;	
	}
	
	if (!cache && _closeishAtoms.size())
	{
		searchAtoms = &_closeishAtoms;
	}

	for (int i = 0; i < searchAtoms->size(); i++)
	{
		AtomPtr search = searchAtoms->at(i);
		if (one == search)
		{
			continue;
		}

		if (!one->closeToAtom(search, tol))
		{
			continue;
		}
		
		std::vector<AtomPtr>::iterator it;
		it = std::find(atomListPtr->begin(), atomListPtr->end(), search);

		if (it != atomListPtr->end())
		{
			continue;	
		}
		
		atomListPtr->push_back(search);
	}
	
	return atoms;
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
		maxTries = 4;
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
	
	addParamType(ParamOptionMaxTries, maxTries);
	maxTries = 1;
	
	while (topAtoms.size() > 0)
	{
		BondPtr topBond;
		int count = 0;

		bool changed = true;

		while (changed && count < maxTries)
		{
			setupNelderMead();
			setCrystal(target);
			setScoreType(scoreType);

			if (rType == RefinementModelPos 
			    || rType == RefinementSavedPos 
			    || rType == RefinementFine 
			    || rType == RefinementCrude 
			    || rType == RefinementMouse)
			{
				setSilent();
			}

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

				bool thorough = true;
				if (rType != RefinementFine && 
				    rType != RefinementMouse &&
				    rType != RefinementCrude)
				{
					thorough = false;
				}

				if (hasParameter(ParamOptionThorough))
				{
					thorough = (getParameter(ParamOptionThorough) > 0);
				}

				topBond = setupThoroughSet(bond, thorough);
			}

			for (size_t l = 0; l < _includeForRefine.size(); l++)
			{
				addSampledAtoms(_includeForRefine[l]);
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
		
		if (hasParameter(ParamOptionTopLevelOnly))
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
	setup_space(&workspace);
	workspace.scoreType = scoreType;
	workspace.crystal = crystal;
	workspace.selectAtoms = shared_from_this();
	workspace.flag = flags;
	workspace.filename = plot;

	return scoreWithMapGeneral(&workspace, plot != "");
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

void AtomGroup::prepareComparisonMap(MapScoreWorkspace *ws, vec3 min, vec3 max)
{
	VagFFTPtr obs = ws->crystal->getFFT();
	mat3x3 m = obs->getRecipBasis();
	mat3x3_mult_vec(m, &min);
	mat3x3_mult_vec(m, &max);
	min -= make_vec3(1, 1, 1);
	max += make_vec3(1, 1, 1);

	VagFFTPtr vag = obs->subFFT(min.x, max.x, min.y, max.y, min.z, max.z);
	ws->comparison = vag;
}

void AtomGroup::addToCubicMap(VagFFTPtr scratchFull, bool saved)
{
	size_t nElements = totalElements();
	scratchFull->prepareAtomSpace();
	_t1->clear();
	_t2->clear();
	_t3->clear();
	
	for (size_t i = 0; i < nElements; i++)
	{
		int count = 0;

		_t1->start();

		for (size_t j = 0; j < atomCount(); j++)
		{
			AtomPtr a = atom(j);
			
			if (a->getElement() != _elements[i])
			{
				continue;
			}

			count++;
			scratchFull->addAtom(a, saved);
		}

		_t1->stop();
	}

	_t3->start();
	scratchFull->fft(FFTAtomsToReal);
	_t3->stop();
	
	for (size_t i = 0; i < _blobs.size(); i++)
	{
		_blobs[i]->addToCubicMap(scratchFull);
	}
}

void AtomGroup::prepareCubicMap(VagFFTPtr *scratchFull, vec3 min, vec3 max, 
                                bool addScratch)
{
	double cubeDim = Options::getActiveCrystal()->getProteinSampling();
	
	if (hasParameter(ParamOptionProteinSampling))
	{
		cubeDim = getParameter(ParamOptionProteinSampling);
	}
	
	/* 2 Angstroms buffer region on either side of the protein */
	double buff = BUFFER_REGION;
	
	if (hasParameter(ParamOptionPadding))
	{
		buff = getParameter(ParamOptionPadding);
	}

	buff += _padding + cubeDim * 3;
	vec3 buffer = make_vec3(buff, buff, buff);
	vec3_subtract_from_vec3(&min, buffer);

	/* Modify the offset to include buffer region */
	vec3 offset = min;
	vec3_add_to_vec3(&max, buffer);
	vec3 limits = vec3_subtract_vec3(max, min);

	/* find out how many voxels in each dimension to cover entire atom group
	 * with a regular grid */
	vec3 extent = limits;
	vec3_mult(&extent, 1 / cubeDim);
	extent.x = floor(extent.x) + 1;
	extent.y = floor(extent.y) + 1;
	extent.z = floor(extent.z) + 1;
	
	CrystalPtr cryst = Options::getActiveCrystal();
	size_t nEle = cryst->totalElements();
	/* two scratches: first will be "constant" FFT and second will
	 * store observed values */
	int scratch = (addScratch ? 1 : 0);
	(*scratchFull) = VagFFTPtr(new VagFFT(extent.x, extent.y, extent.z, 
	                                      nEle, scratch));

	double bFac = cryst->getRealBFactor();

	(*scratchFull)->setBFactor(bFac);
	(*scratchFull)->setScale(cubeDim);
	(*scratchFull)->setOrigin(offset);
	for (int i = 0; i < nEle; i++)
	{
		(*scratchFull)->addElement(cryst->_elements[i]);
	}

	(*scratchFull)->makePlans();
}

void AtomGroup::addToMap(VagFFTPtr fft, bool saved)
{
	size_t nElements = totalElements();
	
	if (nElements == 0)
	{
		return;
	}
	
	VagFFTPtr scratchFull;

	double buffer = BUFFER_REGION;
	vec3 min, max;
	xyzLimits(&min, &max);
	
	for (size_t i = 0; i < _blobs.size(); i++)
	{
		_blobs[i]->getLoopBoundaries(&min, &max);
	}

	prepareCubicMap(&scratchFull, min, max);

	addToCubicMap(scratchFull, saved);

	VagFFT::operation(fft, scratchFull, MapScoreTypeNone, NULL, false);
}

void AtomGroup::createConstantFraction(MapScoreWorkspace *ws,
                                       vec3 min, vec3 max)
{
	/* Now we can add neighbouring constant atoms from the same Crystal
	 * once we have established they do not go over the border. */

	/* Make half-box measures */
	vec3 minmax = vec3_add_vec3(min, max);
	vec3_mult(&minmax, 0.5);
	vec3 half = vec3_subtract_vec3(max, min);
	vec3_mult(&half, 0.5);

	CrystalPtr crystal = ws->crystal;
	ws->extra = crystal->getAtomsInBox(minmax, half.x, half.y, half.z);

	AtomGroupPtr added = AtomGroupPtr(new AtomGroup());

	/* We want to add anything which is static to the 
	 * constant fraction. */
	for (size_t i = 0; i < ws->extra->atomCount(); i++)
	{
		AtomPtr anAtom = ws->extra->atom(i);

		if (ws->selectAtoms->hasAtom(anAtom) ||
		    added->hasAtom(anAtom))
		{
			continue;
		}

		added->addAtom(anAtom);
	}

	ws->extra = added;
}

double AtomGroup::scoreWithReciprocal(MapScoreWorkspace *ws)
{
	if (!ws->selectAtoms->atomCount())
	{
		return 0;
	}

	bool first = (ws->recip == VagFFTPtr());
	
	if (first)
	{
		CrystalPtr crystal = ws->crystal;
		/* first get min/max from entire crystal to make constant
		 * fraction, which is everything we are not */
		vec3 min, max;
		crystal->xyzLimits(&min, &max);
		createConstantFraction(ws, min, max);

		/* copy the exact dimensions of the crystal's FFT */
		ws->recip = VagFFTPtr(new VagFFT(*crystal->getFFT(), 3));
		ws->recip->makePlans();

		/* grab the solvent */
		
		if (crystal->getBucket())
		{
			VagFFTPtr solv = crystal->getBucket()->getSolvent();
			ws->recip->copyFrom(solv);
			ws->recip->copyToScratch(1);
		}
		
		/* get the data amplitudes and copy to scratch */
		DiffractionPtr data = Options::getRuntimeOptions()->getActiveData();
		data->copyToFFT(ws->recip);
		ws->recip->copyToScratch(2);

		/* now we get a new min/max fraction for creation of our
		 * own segment map */
		ws->selectAtoms->xyzLimits(&min, &max);
		ws->selectAtoms->prepareCubicMap(&ws->segment, min, max, false);
	}

	if (first || ws->recalc)
	{
		/* Add the constant atoms to the reciprocal map */
		ws->recip->setStatus(FFTRealSpace);
		ws->recip->multiplyAll(0);
		ws->extra->addToMap(ws->recip);
		ws->recip->fft(FFTRealToReciprocal);

		/* copy the constant fraction into the scratch */
		ws->recip->copyToScratch(0);
		ws->recalc = false;
	}

	double mr = ws->crystal->getMaxResolution();

	ws->recip->multiplyAll(0);
	ws->recip->setStatus(FFTRealSpace);
	ws->selectAtoms->addToCubicMap(ws->segment);

	VagFFT::operation(ws->recip, ws->segment, MapScoreTypeNone);
	ws->recip->fft(FFTRealToReciprocal);
	ws->recip->addScratchBack(0);
	ws->recip->applySymmetry(true, mr);
	ws->recip->addScratchBack(1);

	double cc = ws->recip->compareReciprocalToScratch(2);
	return -cc;
}

double AtomGroup::scoreWithMapGeneral(MapScoreWorkspace *workspace,
                                      bool plot)
{
	workspace->tMap->start();

	if (plot && workspace->filename == "")
	{
		workspace->filename = "cc_score.png";
	}

	if (!workspace->selectAtoms->atomCount())
	{
		return 0;
	}

	AtomGroupPtr selected = workspace->selectAtoms;

	bool first = (workspace->segment == VagFFTPtr());
	
	/* this is the first time we are running the comparison */
	if (first)
	{
		vec3 min, max; 
		selected->xyzLimits(&min, &max);

		/* prepare the size of the maps */
		selected->prepareCubicMap(&workspace->segment, min, max, true);

		/* find all the non-moving atoms */
		createConstantFraction(workspace, min, max);

		/* prepare the selection to be compared */
		selected->prepareComparisonMap(workspace, min, max);
	}
			
	/* in the case where we wish to recalculate constant atom positions */
	if (first || workspace->recalc)
	{
		/* Add the acceptable atoms to the map */
		workspace->extra->addToCubicMap(workspace->segment, true);
		/* copy the constant fraction into the scratch */
		workspace->segment->copyToScratch(0);
		
		workspace->recalc = false;
	}
	
	/* add our moving fraction, updated elsewhere */
	workspace->selectAtoms->addToCubicMap(workspace->segment);
	workspace->tMap->stop();
	
	workspace->tScore->start();
	double score = scoreFinalMap(workspace, plot, first);
	workspace->tScore->stop();

	return score;
}

double AtomGroup::scoreFinalMap(MapScoreWorkspace *ws, bool plot,
                                   bool first)
{
	/* Convert real2Frac to crystal coords to get correct segment
	* of the big real space map. */
	CrystalPtr crystal = ws->crystal;
	double cutoff = MAP_VALUE_CUTOFF;

	std::vector<double> xs, ys, weights;

	VagFFTPtr map = ws->crystal->getFFT();
	bool difference = (ws->flag & MapScoreFlagDifference);
	
	if (difference)
	{
		map = crystal->getDiFFT();
	}
	
	MapScoreType mapType = MapScoreTypeCorrel;
	
	/* Actual variables to compare, correlation calculations */
	ws->segment->copyRealToImaginary();
	/* Add constant fraction to model map */
	ws->segment->addScratchBack(0);
	
	if (first)
	{
		ws->vals.clear();
	}

	int step = 1;
	
	if (ws->selectAtoms->hasParameter(ParamOptionStep))
	{
		step = ws->selectAtoms->getParameter(ParamOptionStep);
		
		if (step <= 0) step = 1;
	}

	VagFFT::operation(ws->comparison, ws->segment, mapType, &ws->vals, 
	                  false, !first, step);

	/* Debugging ... writes cc_score.csv and cc_score.png, csv can be
	* looked at with gnuplot quite nicely.*/
	
//	ws->scoreType = ScoreTypeRFactor;
	double score = 0;
	
	if (ws->scoreType == ScoreTypeCorrel)
	{
		score = -correlation(ws->vals);
	}
	
	if (ws->scoreType == ScoreTypeRFactor)
	{
		score = weighted_r_factor(ws->vals);
	}

	if (plot)
	{
		plotCoordVals(ws->vals, difference, cutoff, ws->filename);
	}

	return score;
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
		
		if (vals[i].weight < cutoff) continue;

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
	double count = 0;
	
	for (int i = 0; i < atomCount(); i++)
	{
		vec3 abs = atom(i)->getAbsolutePosition();
		
		if (!vec3_is_sane(abs))
		{
			continue;
		}
		
		count++;
		vec3_add_to_vec3(&sum, abs);
	}
	
	vec3_mult(&sum, 1 / (double)count);
	
	return sum;
}

vec3 AtomGroup::initialCentroid()
{
	vec3 sum = empty_vec3();
	double count = 0;
	
	for (int i = 0; i < atomCount(); i++)
	{
		if (!atom(i)->isFromPDB())
		{
			continue;
		}

		vec3 abs = atom(i)->getInitialPosition();
		vec3_add_to_vec3(&sum, abs);
		count++;
	}
	
	vec3_mult(&sum, 1/count);
	
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
		
		if (b->hasTwist())
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

void AtomGroup::postParseTidy()
{
	sort();
}

void AtomGroup::sort()
{
	std::sort(_atoms.begin(), _atoms.end(), Atom::greater);
}

void AtomGroup::empty()
{
	_atoms.clear();
}

void AtomGroup::printList()
{
	for (int i = 0; i < atomCount(); i++)
	{
		std::cout << atom(i)->longDesc() << std::endl;
	}
}

bool AtomGroup::isFullyTied()
{
	for (int i = 0; i < atomCount(); i++)
	{
		if (!atom(i)->getModel()->hasExplicitPositions())
		{
			return false;
		}
	}
	
	return true;
}

void AtomGroup::reset()
{
	for (int i = 0; i < atomCount(); i++)
	{
		if (!atom(i)->isSidechain())
		{
			continue;
		}

		if (!atom(i)->getModel()->isBond()) 
		{
			continue;
		}
		
		BondPtr bond = ToBondPtr(atom(i)->getModel());
		
		
		bond->reset();
//		bond->getMajor()->refreshPosition();
	}

	refreshPositions();
}

void AtomGroup::convertWaters()
{
	for (int i = 0; i < atomCount(); i++)
	{
		AtomPtr tmp = atom(i);
		if (tmp->isWater() && tmp->getModel()->isAbsolute())
		{
			SpongePtr nov = SpongePtr(new Sponge(tmp));
			tmp->setModel(nov);
			nov->singleRefine();
			nov->getFinalPositions();
		}
	}
	
	std::cout << "Converted waters." << std::endl;
}

AtomList AtomGroup::atomsSimilarTo(AtomPtr a)
{
	std::string type = a->getAtomName();
	std::string id = a->getChainID();
	int resnum = a->getResidueNum();

	return findAtoms(type, resnum, id);
}
