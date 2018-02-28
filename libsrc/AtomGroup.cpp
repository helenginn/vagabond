//
//  AtomGroup.cpp
//  vagabond
//
//  Created by Helen Ginn on 25/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "AtomGroup.h"
#include "Atom.h"
#include "Element.h"
#include "Bond.h"
#include <sstream>
#include "Crystal.h"
#include <iomanip>
#include "CSV.h"
#include "maths.h"
#include "Shouter.h"
#include "../libccp4/ccp4_spg.h"
#include "FlexRegion.h"
#include "Options.h"

AtomPtr AtomGroup::findAtom(std::string atomType)
{
	for (int i = 0; i < atomCount(); i++)
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

	for (int i = 0; i < atoms.size(); i++)
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

std::map<std::string, int> AtomGroup::conformerMap()
{
	std::map<std::string, int> conformerList;

	for (int i = 0; i < atomCount(); i++)
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
	std::map<std::string, int> conformerList = conformerMap();

	return conformerList.size();
}

std::string AtomGroup::conformer(int i)
{
	if (i > conformerMap().size()) return "";

	std::map<std::string, int> conformerList = conformerMap();
	std::map<std::string, int>::iterator it = conformerList.begin();

	for (int j = 0; j < i; j++) it++;

	return it->first;
}

AtomList AtomGroup::findAtoms(std::string atomType)
{
	AtomList list;

	for (int i = 0; i < atomCount(); i++)
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

	for (int i = 0; i < atomCount(); i++)
	{
		total += atom(i)->getElement()->electronCount();
	}

	return total;
}

std::string AtomGroup::getPDBContribution(PDBType pdbType, CrystalPtr crystal)
{
	std::ostringstream stream;
	int numConf = 0;

	if (!atomCount())
	{
		return "";
	}

	if (pdbType == PDBTypeEnsemble)
	{
		/* Get the total number of conformers to worry about */
		std::vector<BondSample> *samples = atom(0)->getModel()->getManyPositions();

		numConf = samples->size();

		for (int j = 0; j < numConf; j++)
		{
			stream << "MODEL " << std::setw(8) << j << std::setw(66) << " " << std::endl;

			for (int i = 0; i < atomCount(); i++)
			{
				if (!atom(i)->getMonomer())
				{
					continue;
				}

				if (atom(i)->getWeighting() <= 0)
				{
					continue;
				}

				stream << atom(i)->getPDBContribution(j);
			}

			stream << "TER" << std::setw(80) << " " << std::endl;
			stream << "ENDMDL" << std::setw(80) << " " << std::endl;
		}

		return stream.str();
	}

	for (int i = 0; i < atomCount(); i++)
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

void AtomGroup::setUseAbsolute()
{
	for (int i = 0; i < atomCount(); i++)
	{
		atom(i)->setKeepModel();
	}
}

void AtomGroup::addAtomsFrom(AtomGroupPtr child)
{
	for (int i = 0; i < child->atomCount(); i++)
	{
		addAtom(child->atom(i));
	}
}

double AtomGroup::getAverageDisplacement()
{
	double sum = 0;
	double count = 0;

	for (int i = 0; i < atomCount(); i++)
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

	for (int i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getElement()->electronCount() <= 1)
		{
			continue;
		}

		if (initial)
		{
			sum += atom(i)->getInitialBFactor();
			count++;
		}
		else
		{
			if (atom(i)->getModel()->isBond())
			{
				BondPtr bond = ToBondPtr(atom(i)->getModel());
				double val = bond->getMeanSquareDeviation();
				sum += val;
				count++;
			}
		}
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
	for (int i = 0; i < atomCount(); i++)
	{
		atom(i)->getModel()->propagateChange(0);
	}
}

void AtomGroup::refreshPositions(bool quick)
{
	for (int i = 0; i < atomCount(); i++)
	{
		if (!atom(i)) continue;

		if (quick)
		{
			atom(i)->getModel()->propagateChange(0);
			atom(i)->getModel()->getFinalPositions();
			continue;
		}

		atom(i)->getModel()->propagateChange(-1, true);
	}
}


int AtomGroup::totalElectrons(int *fcWeighted)
{
	double sum = 0;
	double weighted = 0;

	for (int i = 0; i < atomCount(); i++)
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
	for (int i = 0; i < atomCount(); i++)
	{
		atom(i)->setWeighting(value);
	}
}

void AtomGroup::resetMagicAxes()
{
}

AtomList AtomGroup::topLevelAtoms()
{
	if (!atomCount()) return AtomList();

	AtomList list;

	for (int i = 0; i < 1; i++)
	{
		std::string conf = conformer(i);
		int j = 0;
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

	for (int i = 0; i < atomCount(); i++)
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
	std::cout << "Refining one residue." << std::endl;
	refine(_target, _rType);
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
		maxTries = 6;
		degrees = 4;
		bondNum = 3;
		refineAngles = false;
		break;

		case RefinementFlexibility:
		scoreType = ScoreTypeModelPos;
		maxTries = 3;
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

	for (int n = 0; n < topAtoms.size(); n++)
	{
		AtomPtr topAtom = topAtoms[n].lock();

		if (n > 0) std::cout << "'" << std::flush;

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

			if (shouldRefineMagicAxis(bond))
			{
				bond->calculateMagicAxis();
			}

			int count = 0;

			BondPtr topBond;

			for (int k = 0; k < bond->downstreamAtomGroupCount(); k++)
			{
				while (rType == RefinementFlexibility && count < maxTries)
				{
					FlexRegion flexer;
					flexer.setup();
					flexer.addBond(ToBondPtr(bond), 8);
					flexer.addSingleBondParameters();
					flexer.sample();
					count++;
				}

				if (rType == RefinementFlexibility)
				{
					rType = RefinementModelPos;
					count = 0;
					maxTries = 60;
				}

				bool changed = true;
				bool addFlex = (rType == RefinementFine);

				while (changed && count < maxTries)
				{
					bond->setActiveGroup(k);
					setupNelderMead();
					setCrystal(target);
					setCycles(16);
					topBond = setupTorsionSet(bond, k, bondNum,
					                          deg2rad(degrees), deg2rad(0.04),
					refineAngles, addFlex);
					setScoreType(scoreType);
					
					for (int l = 0; l < _includeForRefine.size(); l++)
					{
						addSampledAtoms(_includeForRefine[l]);
					}

					if (rType == RefinementModelPos || rType == RefinementFine)
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

	refreshPositions(true);
}

double AtomGroup::scoreWithMap(ScoreType scoreType, CrystalPtr crystal, bool plot)
{
	std::vector<AtomPtr> selected;
	for (int i = 0; i < atomCount(); i++)
	{
		selected.push_back(atom(i));
	}

	return scoreWithMapGeneral(scoreType, crystal, plot, selected);
}

double AtomGroup::scoreWithMapGeneral(ScoreType scoreType, CrystalPtr crystal,
                                      bool plot, std::vector<AtomPtr> selected)
{
	double maxDistance = 0;
	FFTPtr map = crystal->getFFT();
	mat3x3 real2Frac = crystal->getReal2Frac();

	vec3 sum = make_vec3(0, 0, 0);

	/* Find centroid of atom set */
	for (int i = 0; i < selected.size(); i++)
	{
		selected[i]->getModel()->getDistribution(true);
		vec3 offset = selected[i]->getModel()->getAbsolutePosition();
		sum = vec3_add_vec3(sum, offset);	
	}

	vec3_mult(&sum, 1 / (double)selected.size());

	if (sum.x != sum.x)
	{
		return 0;
	}

	vec3 ave = sum;

	/* Find the longest distance from the centroid to determine
	* the max FFT dimensions.*/

	for (int i = 0; i < selected.size(); i++)
	{
		/* Refresh absolute position */
		selected[i]->getModel()->getDistribution();
		vec3 offset = selected[i]->getModel()->getAbsolutePosition();

		vec3 diff = vec3_subtract_vec3(offset, sum);
		double distance = vec3_length(diff);

		if (distance > maxDistance)
		{
			maxDistance = distance;
		}
	}

	/* Calculate appropriate box size and setup FFT */
	double scales = Options::getProteinSampling() / 2;
	int n = 2 * (maxDistance + 3.0) / scales;
	if (n % 2 == 1) n--;

	FFTPtr segment = FFTPtr(new FFT());
	segment->create(n);
	segment->setScales(scales);

	mat3x3 basis = make_mat3x3();
	double toReal = 1 / (scales*(double)n);
	mat3x3_scale(&basis, toReal, toReal, toReal);

	for (int i = 0; i < selected.size(); i++)
	{
		selected[i]->addToMap(segment, basis, ave);
	}

	/* We must fix a problem: because the atoms we add are typically
	* connected to other atoms, the tail-ends of the sampled bit
	* will overlap with other atoms, throwing up high-Fo and low-Fc regions
	* which mess up the correlation. We need to add those extra
	* surrounding atoms, but first we convert all non-zero values to
	* a mask, as we do not want to merely extend the problem to the
	* next atom.*/
	segment->aboveValueToMask(0.1);
	segment->avoidWriteToMaskZero();

	/* Now we can add neighbouring atoms from the same Crystal
	* with impunity, and they won't go over the borderline already
	* established.*/

	std::vector<AtomPtr> extra; // to check for duplicates

	for (int i = 0; i < selected.size(); i++)
	{
		std::vector<AtomPtr> clAtoms = crystal->getCloseAtoms(selected[i], 4.5);

		for (int j = 0; j < clAtoms.size(); j++)
		{
			if (std::find(selected.begin(), selected.end(), clAtoms[j]) !=
			    selected.end())
			{
				continue;
			}

			if (std::find(extra.begin(), extra.end(), clAtoms[j]) !=
			    extra.end())
			{
				continue;
			}

			clAtoms[j]->addToMap(segment, basis, ave);
			extra.push_back(clAtoms[j]);
		}
	}

	/* Convert real2Frac to crystal coords to get correct segment
	* of the big real space map. */
	mat3x3_mult_vec(real2Frac, &ave);

	std::vector<double> xs, ys;
	std::vector<CoordVal> vals;

	double cutoff = FFT::score(map, segment, ave, &vals);

	/* For correlation calculations */
	for (int i = 0; i < vals.size(); i++)
	{
		xs.push_back(vals[i].fo);
		ys.push_back(vals[i].fc);
	}

	double scale = scale_factor_cutoff(xs, ys, cutoff);
	cutoff /= scale;	

	for (int i = 0; i < ys.size(); i++)
	{
		ys[i] /= scale;
		vals[i].fc /= scale;
	}

	/* Debugging ... writes cc_score.csv and cc_score.png, csv can be
	* looked at with gnuplot quite nicely.*/

	if (plot)
	{
		CSVPtr csv = CSVPtr(new CSV(6, "x", "y", "z", "fo", "fc", "mask"));

		for (int i = 0; i < vals.size(); i++)
		{
			double fo = vals[i].fo;
			double fc = vals[i].fc;
			double mask = vals[i].mask;

			vec3 pos = vals[i].pos;
			
			if (fc < cutoff) continue;

			csv->addEntry(6, pos.x, pos.y, pos.z, fo, fc, mask);
		}

		csv->writeToFile("cc_score.csv");

		std::map<std::string, std::string> plotMap;
		plotMap["filename"] = "cc_score";
		plotMap["xHeader0"] = "fc";
		plotMap["yHeader0"] = "fo";
		plotMap["colour0"] = "black";

		plotMap["xTitle0"] = "Fc amplitude";
		plotMap["yTitle0"] = "Fo amplitude";
		plotMap["style0"] = "scatter";
		csv->plotPNG(plotMap);
	}
	
	/* Clear out the massive vectors */
	vals.clear();
	std::vector<CoordVal>().swap(vals);
	/* n.b. this is fucked. please unfuck before continuing. */
	/*    FFT::score(map, obsSeg, zero, NULL, NULL, MapScoreTypeCopyToSmaller);
	*  when you need this ^ check the commits from before 26th Feb. */

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

	return 0;
}

void AtomGroup::addProperties()
{
	addBoolProperty("been_tied", &_beenTied);
	addIntProperty("times_refined", &_timesRefined);

	for (int i = 0; i < atomCount(); i++)
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
