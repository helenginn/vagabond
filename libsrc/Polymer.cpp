//
//  Polymer.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Polymer.h"
#include "Crystal.h"
#include "Sidechain.h"
#include "Monomer.h"
#include "Sampler.h"
#include <iostream>
#include "Backbone.h"
#include "Shouter.h"
#include "Atom.h"
#include "Anchor.h"
#include "Absolute.h"
#include "CSV.h"
#include <sstream>
#include "maths.h"
#include <float.h>
#include "FileReader.h"
#include "Kabsch.h"
#include "Options.h"

void Polymer::addMonomer(MonomerPtr monomer)
{
	if (monomer)
	{
		int resNum = monomer->getResidueNum();
		_monomers[resNum] = monomer;
		monomer->setPolymer(shared_from_this());

		if (monomer->getResidueNum() + 1 > _totalMonomers)
		{
			_totalMonomers = monomer->getResidueNum() + 1;
		}
	}
}

void Polymer::checkChainContinuity()
{
	int foundFirst = -1;
	int foundGap = -1;

	for (int i = 0; i < monomerCount(); i++)
	{
		if (foundFirst >= 0 && foundGap >= 0 && getMonomer(i))
		{
			// Now the monomers have started again

			shout_at_user("Continuity break in chain " + getChainID()
						  + " residues " + i_to_str(foundFirst) + "-"
						  + i_to_str(foundGap) + "\nPlease rebuild and rerun.");
		}

		if (getMonomer(i))
		{
			foundFirst = i;
		}

		if (foundFirst >= 0 && !getMonomer(i))
		{
			foundGap = i;
		}
	}
}

void Polymer::tieAtomsUp()
{
	if (!getMonomer(_anchorNum) || !getMonomer(_anchorNum)->findAtom("N"))
	{
		shout_at_user("Anchor point specified isn't an available residue.\n"\
					  "Please specify an existing residue as an anchor point\n"\
					  "with option --anchor-res=");
	}

	checkChainContinuity();

	AtomPtr n = getMonomer(_anchorNum)->findAtom("N");

	ModelPtr nModel = n->getModel();
	if (nModel->isAbsolute())
	{
		ToAbsolutePtr(nModel)->setBFactor(2.0);
	}

	for (int i = _anchorNum; i < monomerCount(); i++)
	{
		if (getMonomer(i))
		{
			getMonomer(i)->tieAtomsUp();
		}
	}

	for (int i = _anchorNum - 1; i >= 0; i--)
	{
		if (getMonomer(i))
		{
			getMonomer(i)->tieAtomsUp();
		}
	}

	for (int i = 0; i < monomerCount(); i++)
	{
		if (getMonomer(i))
		{
			getMonomer(i)->getSidechain()->splitConformers();

			if (Options::enableTests() && (i == 62 || i == 123))
			{
				std::cout << "Parameterising residue " << i_to_str(i) << std::endl;
				getMonomer(i)->getSidechain()->parameteriseAsRotamers();
			}

			getMonomer(i)->getSidechain()->setInitialDampening();
		}
	}

	resetMagicAxes();
}

void Polymer::summary()
{
	Molecule::summary();
	std::cout << "| I am a polymer with " << _monomers.size() << " monomers." << std::endl;

}

void Polymer::refineMonomer(MonomerPtr monomer, CrystalPtr target,
							RefinementType rType)
{
	if (!monomer)
	{
		return;
	}

	BackbonePtr backbone = monomer->getBackbone();

	if (backbone)
	{
		backbone->refine(target, rType);
	}


	SidechainPtr victim = monomer->getSidechain();

	if (victim && victim->canRefine())
	{
		victim->refine(target, rType);
	}
}

void Polymer::refine(CrystalPtr target, RefinementType rType)
{
	time_t wall_start;
	time(&wall_start);

	std::cout << std::endl;
	std::cout << "Refining chain " << getChainID() << " from anchor to N-terminus...";
	std::cout << std::endl << std::endl << "\t";

	int count = 0;

	for (int i = _anchorNum - 1; i >= 0; i--)
	{
		MonomerPtr monomer = getMonomer(i);

		refineMonomer(monomer, target, rType);

		if (monomer && count % 30 == 30 - 1)
		{
			std::cout << " - " << i << "\n\t" << std::flush;
		}

		count++;
	}

	std::cout << std::endl << std::endl;
	std::cout << "Refining chain " << getChainID() << " from anchor to C-terminus...";
	std::cout << std::endl << std::endl << "\t";

	count = 0;

	for (int i = _anchorNum; i < monomerCount(); i++)
	{
		MonomerPtr monomer = getMonomer(i);
		refineMonomer(monomer, target, rType);

		if (monomer && count % 30 == 30 - 1)
		{
			std::cout << " - " << i << "\n\t" << std::flush;
		}

		count++;
	}

	std::cout << std::endl << std::endl;

	shout_timer(wall_start, "refinement");

}

std::string Polymer::makePDB(PDBType pdbType, CrystalPtr crystal)
{
	std::ostringstream stream;
	stream << getPDBContribution(pdbType, crystal);

	return stream.str();
}

void Polymer::differenceGraphs(std::string graphName, CrystalPtr diffCrystal)
{
	CSVPtr perCA = CSVPtr(new CSV(3, "resnum", "cc", "diffcc"));

	std::vector<double> tempCCs, tempDiffCCs, tempNs;
	double sumCC = 0; double sumDiffCC = 0;

	for (int n = 0; n < monomerCount(); n++)
	{
		if (!getMonomer(n))
		{
			continue;
		}

		AtomPtr ca = getMonomer(n)->getBackbone()->findAtom("CA");
		std::vector<double> xs, ys;
		double cutoff = ca->scoreWithMap(diffCrystal, &xs, &ys, false, MapScoreTypeCorrel);
		double cc = correlation(xs, ys, cutoff);
		sumCC += cc;

		cutoff = ca->scoreWithMap(diffCrystal, &xs, &ys, true, MapScoreTypeCorrel);
		double diffcc = weightedMapScore(xs, ys);
		sumDiffCC += diffcc;

		tempCCs.push_back(cc);
		tempDiffCCs.push_back(diffcc);
		tempNs.push_back(n + 1);
	}

	sumCC /= tempNs.size();
	sumDiffCC /= tempNs.size();

	for (int i = 0; i < tempNs.size(); i++)
	{
		double ccRelative = tempCCs[i];
		double diffCCRelative = tempDiffCCs[i] / sumDiffCC - 1;
		perCA->addEntry(3, tempNs[i], ccRelative, diffCCRelative);
	}

	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = "diffmap_" + graphName;
	plotMap["height"] = "700";
	plotMap["width"] = "1200";
	plotMap["xHeader0"] = "resnum";
	plotMap["yHeader0"] = "cc";
	plotMap["xHeader1"] = "resnum";
	plotMap["yHeader1"] = "diffcc";
	plotMap["yMin0"] = "-1";
	plotMap["yMin1"] = "-1";
	plotMap["yMax0"] = "1";
	plotMap["yMax1"] = "1";
	plotMap["colour0"] = "black";
	plotMap["colour1"] = "red";
	plotMap["xTitle0"] = "Residue number";
	plotMap["yTitle0"] = "Correlation coefficient";
	plotMap["style0"] = "line";
	plotMap["style1"] = "line";

	perCA->writeToFile("diffmap_" + graphName + ".csv");
	perCA->plotPNG(plotMap);
}

void Polymer::graph(std::string graphName)
{
	CSVPtr csv = CSVPtr(new CSV(7, "resnum", "newB", "oldB",
								"pos", "sidepos", "flex", "ellipsoid"));
	CSVPtr csvDamp = CSVPtr(new CSV(4, "resnum", "dN-CA", "dCA-C", "dC-N"));
	CSVPtr csvBlur = CSVPtr(new CSV(4, "resnum", "bN-CA", "bCA-C", "bC-N"));
	CSVPtr sidechainCsv = CSVPtr(new CSV(3, "resnum", "oldB", "newB"));

	for (int i = 0; i < monomerCount(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}

		BackbonePtr backbone = getMonomer(i)->getBackbone();
		SidechainPtr sidechain = getMonomer(i)->getSidechain();
		AtomPtr ca = backbone->findAtom("CA");
		ModelPtr caModel = ca->getModel();
		AtomPtr n = backbone->findAtom("N");
		ModelPtr nModel;
		if (n)
		{
			nModel = n->getModel();
		}

		AtomPtr c = backbone->findAtom("C");
		ModelPtr cModel = c->getModel();

		double value = i;
		double caDampen = 0; double cDampen = 0; double nDampen = 0;
		double caBlur = 0; double cBlur = 0; double nBlur = 0;

		if (caModel->getClassName() == "Bond")
		{
			BondPtr caBond = boost::static_pointer_cast<Bond>(caModel);
			double meanSq = caBond->getMeanSquareDeviation();
			vec3 origEllipsoid = ca->getEllipsoidLongestAxis();
			vec3 newEllipsoid = caBond->longestAxis();
			double angleBetween = vec3_angle_with_vec3(newEllipsoid, origEllipsoid);
			if (angleBetween > deg2rad(90))
			{
				angleBetween -= deg2rad(90);
			}

			double posDisp = ca->posDisplacement();
			double sideDisp = sidechain->getAverageDisplacement();
			double flex = caBond->getFlexibilityPotential();

			csv->addEntry(7, value, meanSq, ca->getInitialBFactor(),
						  posDisp, sideDisp, flex, rad2deg(angleBetween));
			caDampen = Bond::getDampening(&*caBond);
			caBlur = Bond::getTorsionBlur(&*caBond);

			if (caDampen > 0) caBlur = 0;
		}
		else
		{
			csv->addEntry(3, value, ca->getInitialBFactor(), ca->getInitialBFactor());
		}

		if (sidechain)
		{
			double initialBee = sidechain->getAverageBFactor(true);
			double nowBee = sidechain->getAverageBFactor(false);

			sidechainCsv->addEntry(3, value, initialBee, nowBee);
		}

		csvDamp->addEntry(4, value, caDampen, cDampen, nDampen);
		csvBlur->addEntry(4, value, caBlur, cBlur, nBlur);
	}

	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = graphName;
	plotMap["height"] = "700";
	plotMap["width"] = "1200";
	plotMap["xHeader0"] = "resnum";
	plotMap["xHeader1"] = "resnum";
//	plotMap["xHeader2"] = "resnum";
	plotMap["yHeader0"] = "newB";
	plotMap["yHeader1"] = "oldB";
	plotMap["colour0"] = "black";
	plotMap["colour1"] = "red";
	plotMap["yMin0"] = "0";
	plotMap["yMin1"] = "0";
	plotMap["yMax0"] = "40";
	plotMap["yMax1"] = "40";

	plotMap["xTitle0"] = "Residue number";
	plotMap["yTitle0"] = "B factor";
	plotMap["style0"] = "line";
	plotMap["style1"] = "line";

	csv->plotPNG(plotMap);
	csv->writeToFile(graphName + ".csv");

	{
		std::map<std::string, std::string> plotMap;
		plotMap["height"] = "700";
		plotMap["width"] = "1200";
		plotMap["xHeader0"] = "resnum";
		plotMap["colour0"] = "black";
		plotMap["colour1"] = "red";
		plotMap["yMin0"] = "0";
		plotMap["yMax0"] = "180";

		plotMap["xTitle0"] = "Residue number";
		plotMap["style0"] = "line";
		plotMap["filename"] = "ellipsoid_" + graphName;
		plotMap["yHeader0"] = "ellipsoid";
		plotMap["yTitle0"] = "Ellipsoid angle";

		csv->plotPNG(plotMap);
	}

	{
		std::map<std::string, std::string> plotMap;
		plotMap["filename"] = "sidechain_" + graphName;
		plotMap["height"] = "700";
		plotMap["width"] = "1200";
		plotMap["xHeader0"] = "resnum";
		plotMap["yHeader0"] = "newB";
		plotMap["xHeader1"] = "resnum";
		plotMap["yHeader1"] = "oldB";
		plotMap["colour0"] = "black";
		plotMap["colour1"] = "red";
		plotMap["yMin0"] = "0";
		plotMap["yMin1"] = "0";
		plotMap["yMax0"] = "40";
		plotMap["yMax1"] = "40";

		plotMap["xTitle0"] = "Sidechain number";
		plotMap["yTitle0"] = "B factor";
		plotMap["style0"] = "line";
		plotMap["style1"] = "line";

		sidechainCsv->plotPNG(plotMap);
		sidechainCsv->writeToFile("sidechain_" + graphName + ".csv");
	}

	{
		std::map<std::string, std::string> plotMap;
		plotMap["filename"] = "displacement_" + graphName;
		plotMap["height"] = "700";
		plotMap["width"] = "1200";
		plotMap["xHeader0"] = "resnum";
		plotMap["yHeader0"] = "pos";
		plotMap["xHeader1"] = "resnum";
		plotMap["yHeader1"] = "sidepos";
		plotMap["colour0"] = "black";
		plotMap["colour1"] = "blue";
		plotMap["yMin0"] = "0";
		plotMap["yMax0"] = "1";
		plotMap["yMin1"] = "0";
		plotMap["yMax1"] = "1";

		plotMap["xTitle0"] = "Residue number";
		plotMap["yTitle0"] = "Displacement Ang";
		plotMap["style0"] = "line";
		plotMap["style1"] = "line";

		csv->plotPNG(plotMap);
		csv->writeToFile("displacement_" + graphName + ".csv");
	}
}


double Polymer::getSidechainDampening(void *object)
{
	return static_cast<Polymer *>(object)->_sideDampening;
}

void Polymer::setSidechainDampening(void *object, double value)
{
	Polymer *polymer = static_cast<Polymer *>(object);
	polymer->_sideDampening = value;

	for (int i = 0; i < polymer->monomerCount(); i++)
	{
		if (polymer->getMonomer(i))
		{
			polymer->getMonomer(i)->setSidechainDampening(value);
		}
	}
}

double Polymer::getBackboneDampening(void *object)
{
	return static_cast<Polymer *>(object)->_dampening;
}

void Polymer::setBackboneDampening(void *object, double value)
{
	Polymer *polymer = static_cast<Polymer *>(object);
	polymer->_dampening = value;

	for (int i = 0; i < polymer->monomerCount(); i++)
	{
		if (polymer->getMonomer(i))
		{
			polymer->getMonomer(i)->setBackboneDampening(value);
		}
	}
}

void Polymer::findAnchorNearestCentroid()
{
	vec3 sum = make_vec3(0, 0, 0);
	double count = 0;

	for (int i = 0; i < monomerCount(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}

		AtomPtr n = getMonomer(i)->findAtom("N");

		if (!n)
		{
			continue;
		}

		n->getModel()->getFinalPositions();
		vec3 absN = n->getModel()->getAbsolutePosition();

		if (!n) continue;

		count++;
		sum = vec3_add_vec3(sum, absN);
	}

	vec3_mult(&sum, 1 / count);
	int anchorRes = -1;
	double lowestLength = FLT_MAX;

	for (int i = 0; i < monomerCount(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}

		AtomPtr n = getMonomer(i)->findAtom("N");

		if (!n) continue;

		vec3 absCa = n->getModel()->getAbsolutePosition();

		vec3 diff = vec3_subtract_vec3(absCa, sum);
		double length = vec3_length(diff);

		if (length < lowestLength)
		{
			anchorRes = getMonomer(i)->getResidueNum();
			lowestLength = length;
		}
	}

	if (anchorRes < 0)
	{
		shout_at_user("You appear to have no C-alpha atoms in your structure.");
	}
	else
	{
		std::cout << "Anchoring at residue " << anchorRes << ", chain "
		<< getChainID() << std::endl;
	}

	setAnchor(anchorRes);
}

void Polymer::changeAnchor(int num)
{
	resetInitialPositions();

	MonomerPtr newMono = getMonomer(num);

	if (!newMono)
	{
		shout_at_helen("Attempt to anchor residue " + i_to_str(num - 1) +
					   " failed because it doesn't exist on Chain "
					   + getChainID() + ".");
		return;
	}

	AtomPtr newAnchorAtom = newMono->getBackbone()->findAtom("CA");

	if (!newAnchorAtom)
	{
		shout_at_helen("Anchor position CA does not exist\n"\
					   "for residue " + i_to_str(num));
	}

	newMono->getBackbone()->setAnchor();

	for (int i = 0; i < monomerCount(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}

		BackbonePtr bone = getMonomer(i)->getBackbone();
		AtomPtr betaCarbon = bone->betaCarbonTorsionAtom();
		getMonomer(i)->getSidechain()->fixBackboneTorsions(betaCarbon);
		getMonomer(i)->resetMagicAxes();
	}

	_anchorNum = num;
}

void Polymer::setInitialKick(void *object, double value)
{
	Polymer *polymer = static_cast<Polymer *>(object);
	int monomerNum = polymer->_anchorNum;
	polymer->getMonomer(monomerNum)->setKick(value, false);
	polymer->getMonomer(monomerNum - 1)->setKick(value, true);
}

/* For side chains, obviously */
double Polymer::getSideKick(void *object)
{
	return static_cast<Polymer *>(object)->_sideKick;
}

void Polymer::setSideKick(void *object, double value)
{
	Polymer *polymer = static_cast<Polymer *>(object);
	polymer->_sideKick = value;

	for (int i = 0; i < polymer->monomerCount(); i++)
	{
		if (!polymer->getMonomer(i))
		{
			continue;
		}

		polymer->getMonomer(i)->setSideKick(value);
	}
}

double Polymer::getInitialKick(void *object)
{
	Polymer *polymer = static_cast<Polymer *>(object);
	int monomerNum = polymer->_anchorNum;
	return polymer->getMonomer(monomerNum)->getKick();
}

void Polymer::scaleFlexibilityToBFactor(CrystalPtr target)
{
	setupNelderMead();
	setScoreType(ScoreTypeModelFlexiness);
	double value = target->getOverallBFactor();
	std::cout << "Scaling flexibility to " << value << std::endl;
	setTargetFlexibility(value);
//	setVerbose();
	addOverallKickAndDampen(shared_from_this());
	setCycles(50);
	ModelPtr anchor = getAnchorModel();

	addSampledBackbone(shared_from_this());
	setJobName("dampen");
	sample();

	return;
}

ModelPtr Polymer::getAnchorModel()
{
	MonomerPtr anchoredRes = getMonomer(getAnchor());
	ModelPtr model = anchoredRes->findAtom("N")->getModel();

	return model;
}

void Polymer::superimpose()
{
	minimiseRotations();
	minimiseCentroids();
}

void Polymer::minimiseCentroids()
{
	_centroidOffsets.clear();

	std::vector<vec3> addedVecs;
	int count = 0;

	ModelPtr anchor = getAnchorModel();
	getAnchorModel()->getFinalPositions();
	vec3 oldPos = getAnchorModel()->getAbsolutePosition();

	for (int i = 0; i < monomerCount(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}

		count++;

		AtomPtr ca = getMonomer(i)->findAtom("CA");

		if (!ca) continue;

		std::vector<BondSample> someSamples;
		someSamples = ca->getModel()->getFinalPositions();
		std::vector<BondSample> *samples = &someSamples;


		if (!addedVecs.size())
		{
			addedVecs = std::vector<vec3>(samples->size(), make_vec3(0, 0, 0));
		}

		for (int j = 0; j < samples->size(); j++)
		{
			addedVecs[j] = vec3_add_vec3(addedVecs[j], samples->at(j).start);
		}
	}

	// calculate the number of atoms which have been gone into each anchor
	double mult = 1 / (double)count;

	// now we divide every state of the anchor by this.

	for (int i = 0; i < addedVecs.size(); i++)
	{
		vec3_mult(&addedVecs[i], mult);
}

	vec3 meanPos = make_vec3(0, 0, 0); // mean of all centroids

	for (int i = 0; i < addedVecs.size(); i++)
	{
		meanPos = vec3_add_vec3(meanPos, addedVecs[i]);
	}

	vec3_mult(&meanPos, 1 / (double)addedVecs.size());

	// Find the offsets to bring all centroids to the mean value
	for (int i = 0; i < addedVecs.size(); i++)
	{
		vec3 offset = vec3_subtract_vec3(addedVecs[i], meanPos);
		_centroidOffsets.push_back(offset);
	}

	propagateChange();

	anchor->getFinalPositions();
	vec3 newPos = anchor->getAbsolutePosition();

	vec3 backToOld = vec3_subtract_vec3(oldPos, newPos);
	std::cout << "Additional correction: " << vec3_desc(backToOld) << std::endl;

	for (int i = 0; i < _centroidOffsets.size(); i++)
	{
		vec3 offset = vec3_add_vec3(_centroidOffsets[i], backToOld);
		_centroidOffsets[i] = offset;
	}

	std::cout << "Calculated corrections for " << _centroidOffsets.size() << " structures." << std::endl;
}

void Polymer::minimiseRotations()
{
	int num = 0;
	_rotations.clear();
	_centroids.clear();

	for (int i = 0; i < monomerCount(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}

		AtomPtr ca = getMonomer(i)->findAtom("CA");
		if (!ca) continue;

		std::vector<BondSample> *samples;
		samples = ca->getModel()->getManyPositions(BondSampleThorough);

		num = samples->size();

		break;
	}

	int f = num / 2;
	std::vector<mat3x3> tmpMats;

	for (int i = 0; i < num; i++)
	{
		std::vector<double> weights;
		std::vector<vec3> fixedVecs, variantVecs;

		for (int j = 0; j < monomerCount(); j++)
		{
			if (!getMonomer(j))
			{
				continue;
			}

			BackbonePtr bone = getMonomer(j)->getBackbone();

			for (int k = 0; k < bone->atomCount(); k++)
			{
				AtomPtr anAtom = bone->atom(k);
				if (!anAtom) continue;

				std::vector<BondSample> *samples;
				ModelPtr model = anAtom->getModel();

				if (!model)
				{
					shout_at_helen("Missing model for backbone atom!");
				}

				double bee = model->getMeanSquareDeviation();
				double weight = 1 / (bee * bee);
				weights.push_back(weight);

				samples = model->getManyPositions(BondSampleThorough);
				vec3 fixed = model->getAbsolutePosition();

				vec3 variant = samples->at(i).start;
				fixedVecs.push_back(fixed);
				variantVecs.push_back(variant);
			}
		}

		Kabsch kabsch;

		kabsch.setAtoms(fixedVecs, variantVecs);
		kabsch.setWeights(weights);
		_centroids.push_back(kabsch.fixCentroids());
		mat3x3 mat = kabsch.run();

		if (kabsch.didFail())
		{
			_centroidOffsets.clear();
			return;
		}

		tmpMats.push_back(mat);
	}

	_rotations = tmpMats;

	std::cout << "Kabsch'd them all!" << std::endl;

	propagateChange();
}

void Polymer::closenessSummary()
{
	double posSum = 0;
	double bSum = 0;
	int count = 0;

	for (int i = 0; i < atomCount(); i++)
	{
		double disp = atom(i)->posDisplacement();

		if (disp != disp) continue;

		posSum += disp;
		bSum += atom(i)->getModel()->getMeanSquareDeviation();
		count++;
	}

	posSum /= count;
	bSum /= count;

	std::cout << "Across all Chain " << getChainID() << " atoms:\n";
	std::cout << "\tB factor (Å^2): " << bSum << std::endl;
	std::cout << "\tPositional displacement from PDB (Å): " << posSum << std::endl;
}

void Polymer::downWeightResidues(int start, int end, double value) // inclusive
{
	double count = 0;
	for (int i = start; i <= end; i++)
	{
		if (getMonomer(i))
		{
			getMonomer(i)->setWeighting(value);
			count++;
		}
	}

	std::cout << "Set " << count << " residues in region " << getChainID()
	<< start << "-" << end << " to weighting of " << 0 << std::endl;
}

bool Polymer::test()
{
	bool bondsOk = true;

	for (int i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getModel()->isBond())
		{
			bondsOk *= ToBondPtr(atom(i)->getModel())->test();
		}
	}

	if (!bondsOk)
	{
		std::cout << "Bonds FAILED test for polymer " << getChainID() << std::endl;
	}

	return true;
}

void Polymer::reportParameters()
{
	int count = 0;

	for (int i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getModel()->isBond())
		{
			if (ToBondPtr(atom(i)->getModel())->isRefinable())
			{
				count++;
			}
		}
	}

	std::cout << "Chain " << getChainID() << " has " << count
	<< " refinable bonds." << std::endl;
}
