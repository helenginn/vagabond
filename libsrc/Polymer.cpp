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
#include "FlexGlobal.h"
#include "RefinementNelderMead.h"

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
		ToAbsolutePtr(nModel)->setBFactor(_startB);
		ToAbsolutePtr(nModel)->setAnchorPoint();
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

	double kick = Options::getKick();
	setInitialKick(this, kick);

	for (int i = 0; i < monomerCount(); i++)
	{
		if (getMonomer(i))
		{
			getMonomer(i)->getSidechain()->splitConformers();

			if ((Options::enableTests() == 1 || Options::enableTests() == 2)
					&& (i == 62 || i == 30 || i == 78 ||
						i >= 123 || i == 103))
			{
				getMonomer(i)->getSidechain()->parameteriseAsRotamers();
			}
			/*
			   if ((Options::enableTests() == 3)
			   && (i == 857 || i == 962 || i == 1011 || i == 1012 ||
			   i == 1015 || i == 1028 || i == 1053 || i == 1068))
			   {
			   getMonomer(i)->getSidechain()->parameteriseAsRotamers();
			   }
			   */

			getMonomer(i)->getSidechain()->setInitialDampening();
		}
	}

	resetMagicAxes();
}

void Polymer::splitConformers()
{
	for (int i = 0; i < monomerCount(); i++)
	{
		if (getMonomer(i))
		{
			getMonomer(i)->getSidechain()->splitConformers();
		}
	}
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

void Polymer::refineToEnd(int monNum, CrystalPtr target, RefinementType rType)
{
	int start = monNum;
	int end = (monNum < _anchorNum) ? -1 : monomerCount();
	int skip = (monNum < _anchorNum) ? -1 : 1;

	int count = 0;
	double startCCAve = 0;
	std::map<MonomerPtr, double> preScores;
	
	std::cout << "Refining chain " << getChainID() << " from anchor to ";
	std::cout << (skip > 0 ? "C" : "N");
	std::cout <<  "-terminus...";
	std::cout << "\t";

	for (int i = start; i != end; i += skip)
	{
		MonomerPtr monomer = getMonomer(i);
		if (!monomer)
		{
			continue;
		}
		
		double score = monomer->scoreWithMap(ScoreTypeCorrel, target);
		startCCAve += score;

		preScores[monomer] = score;
		count++;
	}

	startCCAve /= (double)count;
	count = 0;
	
	double endCCAve = 0;

	for (int i = start; i != end; i += skip)
	{
		MonomerPtr monomer = getMonomer(i);
		if (!monomer)
		{
			continue;
		}

		copyParams(monomer->getSidechain());
		copyParams(monomer->getBackbone());
		refineMonomer(monomer, target, rType);

		double score = monomer->scoreWithMap(ScoreTypeCorrel, target);	
		double pre = preScores[monomer];
		
		/* Scores are negative by default */
		double diff = (score - pre) * 100.;
		/* improvement of 2% will be 20 + symbols after the residue! */
		int signs = fabs(diff * 10);
		int dir = (diff < 0);	
		
		std::cout << " ";
		for (int j = 0; j < signs; j++)
		{
			std::cout << (dir ? "+" : "-");
		}
		std::cout << std::endl << "\t";
		
		endCCAve += score;
		count++;
	}

	endCCAve /= (double)count;
	
	std::cout << "Average CC of monomers went ";
	std::cout << ((endCCAve < startCCAve) ? "up " : "down ");
	std::cout << "from " << -startCCAve * 100 << " to " << -endCCAve * 100;
	std::cout << "." << std::endl;

	clearParams();
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
	return;
	CSVPtr perCA = CSVPtr(new CSV(3, "resnum", "cc", "diffcc"));

	std::vector<double> tempCCs, tempDiffCCs, tempNs;
	double sumCC = 0; double sumDiffCC = 0;
	FFTPtr fft = diffCrystal->getFFT();
	FFTPtr difft = diffCrystal->getDiFFT();
	mat3x3 real2Frac = diffCrystal->getReal2Frac();

	for (int n = 0; n < monomerCount(); n++)
	{
		if (!getMonomer(n) || !getMonomer(n)->getBackbone())
		{
			continue;
		}

		BackbonePtr backbone = getMonomer(n)->getBackbone();
		double cc = -getMonomer(n)->scoreWithMap(ScoreTypeCorrel, diffCrystal);
		sumCC += cc;

		double diffcc = -getMonomer(n)->scoreWithMap(ScoreTypeMultiply,
				diffCrystal);
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
	plotMap["yMin0"] = "0";
	plotMap["yMax0"] = "1";
	plotMap["yMin1"] = "-10";
	plotMap["yMax1"] = "10";

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
	//    plotMap["xHeader2"] = "resnum";
	plotMap["yHeader0"] = "newB";
	plotMap["yHeader1"] = "oldB";
	plotMap["colour0"] = "black";
	plotMap["colour1"] = "red";
	plotMap["yMin0"] = "0";
	plotMap["yMin1"] = "0";
	plotMap["yMax0"] = "60";
	plotMap["yMax1"] = "60";

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
	int anchorNum = polymer->_anchorNum;

	for (int i = 0; i < polymer->monomerCount(); i++)
	{
		MonomerPtr monomer = polymer->getMonomer(i);

		if (!monomer) continue;

		BackbonePtr bone = monomer->getBackbone();
		double kick = value * (monomer->getResidueNum() < anchorNum ? -1 : 1);

		for (int j = 0; j < bone->atomCount(); j++)
		{
			ModelPtr model = bone->atom(j)->getModel();

			if (model->isBond())
			{
				BondPtr bond = ToBondPtr(model);
				Bond::setTorsionBlur(&*bond, kick);
			}
		}
	}
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

ModelPtr Polymer::getAnchorModel()
{
	MonomerPtr anchoredRes = getMonomer(getAnchor());
	if (!anchoredRes)
	{
		return ModelPtr();
	}

	ModelPtr model = anchoredRes->findAtom("N")->getModel();

	return model;
}

std::vector<vec3> Polymer::getAnchorSphereDiffs()
{
	std::vector<vec3> results;
	ModelPtr anchor = getAnchorModel();

	if (!anchor)
	{
		return std::vector<vec3>();
	}

	std::vector<BondSample> *finals = getAnchorModel()->getManyPositions();
	vec3 sum = make_vec3(0, 0, 0);

	for (int i = 0; i < finals->size(); i++)
	{
		sum = vec3_add_vec3(sum, finals->at(i).start);
	}

	vec3_mult(&sum, 1/(double)finals->size());

	for (int i = 0; i < finals->size(); i++)
	{
		vec3 onePos = finals->at(i).start;
		vec3 diff = vec3_subtract_vec3(onePos, sum);
		results.push_back(diff);
	}

	return results;
}

void Polymer::applyTranslationTensor()
{
	_transTensorOffsets.clear();
	_extraRotationMats.clear();

	std::vector<vec3> sphereDiffs = getAnchorSphereDiffs();

	for (int i = 0; i < sphereDiffs.size(); i++)
	{
		vec3 diff = sphereDiffs[i];
		vec3 diffTensored = diff;
		mat3x3_mult_vec(_transTensor, &diffTensored);
		vec3 movement = vec3_subtract_vec3(diffTensored, diff);
		_transTensorOffsets.push_back(movement);
	}
}

void Polymer::calculateExtraRotations()
{
	// We have our sphereDiffs.
	std::vector<vec3> sphereDiffs = getAnchorSphereDiffs();

	// We have a magic axis, _magicRotAxis. We want to know
	// how close each difference is to one of the poles, mit sign.

	if (!sphereDiffs.size())
	{
		return;
	}

	for (int i = 0; i < sphereDiffs.size(); i++)
	{
		vec3 diff = sphereDiffs[i];
		double cosine = vec3_cosine_with_vec3(diff, _magicRotAxis);
		vec3 cross = vec3_cross_vec3(diff, _magicRotAxis);

		// Need to make sure the poles are opposite. 
		// Choice of doing so is relatively arbitrary.
		double angleMult = cosine;
		if (cross.x < 0)
		{
			angleMult *= -1;
		}

		// Now we make a matrix around our other vector,
		// _rotationAxis, with the modulated _rotationAngle         

		double angle = _rotationAngle * angleMult;
		mat3x3 rot = mat3x3_unit_vec_rotation(_rotationAxis, angle);
		_extraRotationMats.push_back(rot);
	}    
}

void Polymer::superimpose()
{
	_centroids.clear();
	_centroidOffsets.clear();
	_rotations.clear();
	_transTensorOffsets.clear();
	_extraRotationMats.clear();
	propagateChange();

	minimiseRotations();
	minimiseCentroids();

	ModelPtr model = getAnchorModel();
	AtomPtr ca = getMonomer(monomerCount() - 1)->findAtom("CA");
	ModelPtr one = ca->getModel();

	if (ca)
	{
		AtomPtr nz = getMonomer(monomerCount() - 1)->findAtom("NZ");

		if (nz) one = nz->getModel();

		CSVPtr csv = CSVPtr(new CSV(3, "x", "y", "z"));
		std::vector<vec3> poss;
		poss = one->polymerCorrectedPositions();

		for (int i = 0; i < poss.size(); i++)
		{
			vec3 pos = poss[i];
			csv->addEntry(3, pos.x, pos.y, pos.z);
		}

		csv->writeToFile("ca_pos.csv");
	}

	if (model->isAbsolute())
	{
		std::vector<vec3> sphereAngles = ToAbsolutePtr(model)->getSphereAngles();
		CSVPtr csv = CSVPtr(new CSV(10, "psi", "phi", "theta", "corr_x", "corr_y", "corr_z", "rot_angle", "rot_axis_x", "rot_axis_y", "rot_axis_z"));
		CSVPtr one = CSVPtr(new CSV(3, "x", "y", "z"));

		for (int i = 0; i < sphereAngles.size(); i++)
		{
			mat3x3 rotMat = getRotationCorrections()[i];
			double angle = mat3x3_rotation_angle(rotMat);
			vec3 axis = mat3x3_rotation_axis(rotMat);

			double xMove = getCentroidOffsets()[i].x;
			double yMove = getCentroidOffsets()[i].y;
			double zMove = getCentroidOffsets()[i].z;
			csv->addEntry(10, sphereAngles[i].x, sphereAngles[i].y,
					sphereAngles[i].z, xMove, yMove, zMove, angle,
					axis.x, axis.y, axis.z); 
		}

		csv->writeToFile("kabsch.csv");
	}

	applyTranslationTensor();
}

void Polymer::minimiseCentroids()
{
	std::vector<vec3> addedVecs;
	int count = 0;

	ModelPtr anchor = getAnchorModel();

	if (!anchor)
	{
		return;
	}

	getAnchorModel()->getFinalPositions();
	vec3 oldPos = getAnchorModel()->getAbsolutePosition();

	for (int i = 0; i < monomerCount(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}

		AtomPtr ca = getMonomer(i)->findAtom("CA");

		if (!ca) continue;

		count++;

		std::vector<BondSample> *samples;
		std::vector<BondSample> someSamples;
		someSamples = ca->getModel()->getFinalPositions();
		samples = &someSamples;

		if (!addedVecs.size())
		{
			addedVecs = std::vector<vec3>(samples->size(), make_vec3(0, 0, 0));;
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

	// starting point for rotation centre of whole molecule movements
	// Used by Model to correct for translation stuff

	vec3_mult(&meanPos, 1 / (double)addedVecs.size());
	_rotationCentre = meanPos;

	// Remove the old centroid positions if not already
	_centroidOffsets.clear();

	// Find the offsets to bring all centroids to the mean value
	for (int i = 0; i < addedVecs.size(); i++)
	{
		vec3 offset = vec3_subtract_vec3(addedVecs[i], meanPos);
		_centroidOffsets.push_back(offset);
	}

	propagateChange();

	// This time, it should apply the offset to the anchor.
	// But we want the anchor to remain in the same place as before.
	anchor->getFinalPositions();
	vec3 newPos = anchor->getAbsolutePosition();

	// Find the difference.
	vec3 backToOld = vec3_subtract_vec3(oldPos, newPos);

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

	std::vector<vec3> tmpCentroids;

	for (int i = 0; i < monomerCount(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}

		AtomPtr ca = getMonomer(i)->findAtom("CA");
		if (!ca) continue;

		std::vector<BondSample> *samples;
		samples = ca->getModel()->getManyPositions();

		num = samples->size();

		break;
	}

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

				samples = model->getManyPositions();
				vec3 fixed = model->getAbsolutePosition();

				vec3 variant = samples->at(i).start;
				fixedVecs.push_back(fixed);
				variantVecs.push_back(variant);
			}
		}

		Kabsch kabsch;

		kabsch.setAtoms(fixedVecs, variantVecs);
		kabsch.setWeights(weights);
		tmpCentroids.push_back(kabsch.fixCentroids());
		mat3x3 mat = kabsch.run();

		if (kabsch.didFail())
		{
			_centroidOffsets.clear();
			return;
		}

		tmpMats.push_back(mat);
	}

	_centroids = tmpCentroids;
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
		if (!atom(i)->isFromPDB())
		{
			continue;
		}

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

	std::cout << "Testing polymer... silence is good." << std::endl;

	for (int i = 0; i < atomCount(); i++)
	{
		if (!atom(i)->getModel())
		{
			std::cout << "Missing model for " << atom(i)->shortDesc() << std::endl;
			bondsOk = 0;
		}

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

void Polymer::optimiseWholeMolecule(bool translation, bool rotation)
{
	std::cout << "Optimising whole molecule shifts to match the PDB file." << std::endl;

	NelderMeadPtr nelderMead = NelderMeadPtr(new NelderMead());

	if (translation)
	{
		nelderMead->addParameter(this, getTransTensor11, setTransTensor11, 0.5, 0.01);
		nelderMead->addParameter(this, getTransTensor12, setTransTensor12, 0.1, 0.01);
		nelderMead->addParameter(this, getTransTensor21, setTransTensor21, 0.1, 0.01);
		nelderMead->addParameter(this, getTransTensor13, setTransTensor13, 0.1, 0.01);
		nelderMead->addParameter(this, getTransTensor31, setTransTensor31, 0.1, 0.01);
		nelderMead->addParameter(this, getTransTensor22, setTransTensor22, 0.5, 0.01);
		nelderMead->addParameter(this, getTransTensor23, setTransTensor23, 0.1, 0.01);
		nelderMead->addParameter(this, getTransTensor32, setTransTensor32, 0.1, 0.01);
		nelderMead->addParameter(this, getTransTensor33, setTransTensor33, 0.5, 0.01);
	}

	if (rotation)
	{
		nelderMead->addParameter(this, getRotPhi, setRotPhi, 0.1, 0.0001);
		nelderMead->addParameter(this, getRotPsi, setRotPsi, 0.1, 0.0001);
		nelderMead->addParameter(this, getRotAngle, setRotAngle, 0.01, 0.0001);
		nelderMead->addParameter(this, getRotCentreX, setRotCentreX, 5.0, 0.01);
		nelderMead->addParameter(this, getRotCentreY, setRotCentreY, 5.0, 0.01);
		nelderMead->addParameter(this, getRotCentreZ, setRotCentreZ, 5.0, 0.01);
	}

	nelderMead->setCycles(30);
	nelderMead->setVerbose(true);

	double bFac = getAverageBFactor();

	FlexGlobal target;
	target.setAtomGroup(AtomGroup::shared_from_this());
	target.matchOriginalBees();
	nelderMead->setEvaluationFunction(FlexGlobal::score, &target);
	nelderMead->refine();
}

void Polymer::addProperties()
{
	Molecule::addProperties();

	addIntProperty("anchor_res", &_anchorNum);
	addMat3x3Property("trans_tensor", &_transTensor);

	for (int i = 0; i < monomerCount(); i++)
	{
		if (!getMonomer(i)) continue;

		addChild("monomer", getMonomer(i));
	}
}

void Polymer::addObject(ParserPtr object, std::string category)
{
	if (category == "monomer")
	{
		MonomerPtr monomer = ToMonomerPtr(object);
		addMonomer(monomer);
	}

	Molecule::addObject(object, category);
}

void Polymer::postParseTidy()
{
	Molecule::postParseTidy();
	applyTranslationTensor();
}
