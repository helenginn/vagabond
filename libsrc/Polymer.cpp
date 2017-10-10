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
#include "Bond.h"
#include "Absolute.h"
#include "CSV.h"
#include <fstream>
#include "maths.h"
#include <float.h>
#include "FileReader.h"

void Polymer::addMonomer(MonomerPtr monomer)
{
	long existingMonomers = monomerCount();
	_monomers[existingMonomers] = monomer;
	if (monomer)
	{
		monomer->setPolymer(shared_from_this());
	}
}

void Polymer::tieAtomsUp()
{
	if (!getMonomer(_anchorNum - 1))
	{
		shout_at_user("Anchor point specified isn't an available residue.\n"\
					  "Please specify an existing residue as an anchor point\n"\
					  "with option --anchor-res=");
	}

	for (int i = _anchorNum - 1; i < monomerCount(); i++)
	{
		if (getMonomer(i))
		{
			getMonomer(i)->tieAtomsUp();
		}
	}

	for (int i = _anchorNum - 2; i >= 0; i--)
	{
		if (getMonomer(i))
		{
			getMonomer(i)->tieAtomsUp();
		}
	}
}

void Polymer::summary()
{
	Molecule::summary();
	std::cout << "| I am a polymer with " << monomerCount() << " monomers." << std::endl;

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

	for (int i = _anchorNum - 2; i >= 0; i--)
	{
		MonomerPtr monomer = getMonomer(i);
		refineMonomer(monomer, target, rType);
	}

	for (int i = _anchorNum - 1; i < monomerCount(); i++)
	{
		MonomerPtr monomer = getMonomer(i);
		refineMonomer(monomer, target, rType);
	}

	shout_timer(wall_start, "refinement");

}

void Polymer::makePDB(std::string filename)
{
	std::ofstream file;
	file.open(filename.c_str());

	for (int i = 0; i < monomerCount(); i++)
	{
		MonomerPtr monomer = getMonomer(i);

		if (!monomer)
		{
			continue;
		}

		SidechainPtr victim = monomer->getSidechain();

		file << monomer->getBackbone()->getPDBContribution();
		file << victim->getPDBContribution();
	}

	file.close();

	std::cout << "Written PDB to " << filename << "." << std::endl;
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

	std::cout << "Average CC " << sumCC << std::endl;

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

	std::cout << "Written out " << graphName << " graph set." << std::endl;
}

void Polymer::graph(std::string graphName)
{
	CSVPtr csv = CSVPtr(new CSV(10, "resnum", "newB", "oldB", "oldX", "oldY", "oldZ",
								"newX", "newY", "newZ", "pos"));
	CSVPtr csvDamp = CSVPtr(new CSV(4, "resnum", "dN-CA", "dCA-C", "dC-N"));
	CSVPtr csvBlur = CSVPtr(new CSV(4, "resnum", "bN-CA", "bCA-C", "bC-N"));
	CSVPtr sidechainCsv = CSVPtr(new CSV(3, "resnum", "oldB", "newB", "pos"));

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
		ModelPtr nModel = n->getModel();
		AtomPtr c = backbone->findAtom("C");
		ModelPtr cModel = c->getModel();

		double value = i;
		double caDampen = 0; double cDampen = 0; double nDampen = 0;
		double caBlur = 0; double cBlur = 0; double nBlur = 0;

		if (caModel->getClassName() == "Bond")
		{
			BondPtr caBond = std::static_pointer_cast<Bond>(caModel);
			double meanSq = caBond->getMeanSquareDeviation();
			double newX = caBond->getMeanSquareDeviation(-1, 0);
			double newY = caBond->getMeanSquareDeviation(-1, 1);
			double newZ = caBond->getMeanSquareDeviation(-1, 2);
			double oldX = ca->getInitialAnisoB(0);
			double oldY = ca->getInitialAnisoB(1);
			double oldZ = ca->getInitialAnisoB(2);
			double posDisp = ca->posDisplacement();

			csv->addEntry(10, value, meanSq, ca->getInitialBFactor(),
						  oldX, oldY, oldZ, newX, newY, newZ, posDisp);
			caDampen = Bond::getDampening(&*caBond);
			caBlur = Bond::getTorsionBlur(&*caBond);

			if (caDampen > 0) caBlur = 0;
		}
		else
		{
			double oldX = ca->getInitialAnisoB(0);
			double oldY = ca->getInitialAnisoB(1);
			double oldZ = ca->getInitialAnisoB(2);

			csv->addEntry(9, value, ca->getInitialBFactor(), ca->getInitialBFactor(),
						  oldX, oldY, oldZ, oldX, oldY, oldZ);
		}

		if (cModel->getClassName() == "Bond")
		{
			BondPtr cBond = std::static_pointer_cast<Bond>(cModel);
			cDampen = Bond::getDampening(&*cBond);
			cBlur = Bond::getTorsionBlur(&*cBond);

			if (cDampen > 0) cBlur = 0;
		}

		if (nModel->getClassName() == "Bond")
		{
			BondPtr nBond = std::static_pointer_cast<Bond>(nModel);
			nDampen = Bond::getDampening(&*nBond);
			nBlur = Bond::getTorsionBlur(&*nBond);
			if (nDampen > 0) nBlur = 0;
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
	plotMap["yHeader0"] = "newB";
	plotMap["xHeader1"] = "resnum";
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
		plotMap["xHeader1"] = "resnum";
		plotMap["colour0"] = "black";
		plotMap["colour1"] = "red";
		plotMap["yMin0"] = "0";
		plotMap["yMin1"] = "0";
		plotMap["yMax0"] = "40";
		plotMap["yMax1"] = "40";

		plotMap["xTitle0"] = "Residue number";
		plotMap["style0"] = "line";
		plotMap["style1"] = "line";


		plotMap["filename"] = "aniso_x_" + graphName;
		plotMap["yHeader0"] = "newX";
		plotMap["yHeader1"] = "oldX";
		plotMap["yTitle0"] = "Aniso Bx factor";

		csv->plotPNG(plotMap);

		plotMap["filename"] = "aniso_y_" + graphName;
		plotMap["yHeader0"] = "newY";
		plotMap["yHeader1"] = "oldY";
		plotMap["yTitle0"] = "Aniso By factor";

		csv->plotPNG(plotMap);

		plotMap["filename"] = "aniso_z_" + graphName;
		plotMap["yHeader0"] = "newZ";
		plotMap["yHeader1"] = "oldZ";
		plotMap["yTitle0"] = "Aniso Bz factor";

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
		plotMap["colour0"] = "black";
		plotMap["yMin0"] = "0";
		plotMap["yMax0"] = "1";

		plotMap["xTitle0"] = "Residue number";
		plotMap["yTitle0"] = "Displacement Ang";
		plotMap["style0"] = "line";

		csv->plotPNG(plotMap);
		csv->writeToFile("displacement_" + graphName + ".csv");
	}

	plotMap["filename"] = "dampening_" + graphName;
	plotMap["yHeader0"] = "dN-CA";
	plotMap["xHeader1"] = "resnum";
	plotMap["yHeader1"] = "dCA-C";
	plotMap["xHeader2"] = "resnum";
	plotMap["yHeader2"] = "dC-N";
	plotMap["style2"] = "line";
	plotMap["colour2"] = "blue";
	plotMap["yTitle0"] = "Dampening factor";
	plotMap["yMin0"] = "-2";
	plotMap["yMin1"] = "-2";
	plotMap["yMin2"] = "-2";
	plotMap["yMax0"] = "2";
	plotMap["yMax1"] = "2";
	plotMap["yMax2"] = "2";

	csvDamp->plotPNG(plotMap);

	plotMap["filename"] = "blurring_" + graphName;
	plotMap["yHeader0"] = "bN-CA";
	plotMap["yHeader1"] = "bCA-C";
	plotMap["yHeader2"] = "bC-N";
	plotMap["yTitle0"] = "Blurring factor";
	plotMap["yMin0"] = "-0.5";
	plotMap["yMin1"] = "-0.5";
	plotMap["yMin2"] = "-0.5";
	plotMap["yMax0"] = "0.5";
	plotMap["yMax1"] = "0.5";
	plotMap["yMax2"] = "0.5";

	csvBlur->plotPNG(plotMap);

	std::cout << "Written out " << graphName << " graph set." << std::endl;
}

double Polymer::getConstantDampening(void *object)
{
	return static_cast<Polymer *>(object)->_dampening;
}

void Polymer::setConstantDampening(void *object, double value)
{
	Polymer *polymer = static_cast<Polymer *>(object);
	polymer->_dampening = value;

	for (int i = 0; i < polymer->monomerCount(); i++)
	{
		if (polymer->getMonomer(i))
		{
			polymer->getMonomer(i)->setConstantDampening(value);
		}
	}
}

void Polymer::changeAnchor(int num)
{
	resetInitialPositions();

	int oldAnchor = _anchorNum;

	MonomerPtr newMono = getMonomer(num - 1);

	if (!newMono)
	{
		shout_at_helen("Monomer " + i_to_str(num - 1) + " doesn't exist.");
		return;
	}

	AtomPtr newAnchorAtom = newMono->getBackbone()->findAtom("CA");

	if (!newAnchorAtom)
	{
		shout_at_helen("Anchor position CA does not exist\n"\
					   "for residue " + i_to_str(num - 1));
	}

	/* Tasks:
	 * - convert old anchor to a normal bond
	 * (should be in direction of tying - polarity of (oldAnchor - newAnchor)
	 * - convert old bond to a new anchor
	 * - write the reversal
	 */

	newMono->getBackbone()->setAnchor();

	int limit = (oldAnchor < num) ? 0 : monomerCount();
	int step = (oldAnchor < num) ? -1 : 1;

	for (int i = num - 1; i < limit; i += step)
	{
		if (!getMonomer(i))
		{
			return;
		}

		BackbonePtr bone = getMonomer(i)->getBackbone();
		AtomPtr betaCarbon = bone->betaCarbonTorsionAtom();
		getMonomer(i)->getSidechain()->fixBackboneTorsions(betaCarbon);
	}

	_anchorNum = num;
}

void Polymer::setInitialKick(void *object, double value)
{
	Polymer *polymer = static_cast<Polymer *>(object);
	int monomerNum = polymer->_anchorNum - 1;
	polymer->getMonomer(monomerNum)->setKick(value, 0);
	polymer->getMonomer(monomerNum - 1)->setKick(value, 1);
}

double Polymer::getInitialKick(void *object)
{
	Polymer *polymer = static_cast<Polymer *>(object);
	int monomerNum = polymer->_anchorNum - 1;
	return polymer->getMonomer(monomerNum)->getKick();
}

void Polymer::scaleFlexibilityToBFactor(CrystalPtr target)
{
	setupNelderMead();
	setScoreType(ScoreTypeModelOverallB);
	double value = target->getOverallBFactor();
	setOverallBFactor(value);
	addOverallKickAndDampen(shared_from_this());
	setCycles(50);

	AtomPtr nitrogen = getMonomer(getAnchor() - 1)->findAtom("N");

	if (!nitrogen)
	{
		shout_at_helen("Mistake");
	}

//	addAbsoluteBFactor(ToAbsolutePtr(nitrogen->getModel()), 0.4, 0.1);
	addSampledBackbone(shared_from_this());
	setJobName("kick_and_dampen");
	sample();
}
