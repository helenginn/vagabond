//
//  Polymer.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Polymer.h"
#include "BoneDensity.h"
#include "Timer.h"
#include "Crystal.h"
#include "Sidechain.h"
#include "Monomer.h"
#include "Sampler.h"
#include <iostream>
#include "Backbone.h"
#include "Shouter.h"
#include "Atom.h"
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
#include "RefinementGridSearch.h"
#include "Hydrogenator.h"

Polymer::Polymer()
{
	_dampening = Options::getDampen();
	_kick = Options::getKick();
	_sideDampening = 0.05;
	_sideKick = 0;
	_anchorNum = 0;
	_totalMonomers = 0;
	_transTensor = make_mat3x3();
	_transExponent = 0;
	_rotExponent = 0;
	mat3x3_scale(&_transTensor, 1.5, 1.5, 1.5);
	_overallScale = 0;
	_startB = Options::getBStart();
	_extraRotParams = {1, 0, 0};
	_tmpPhi = 0;
	_tmpPsi = 0;
}

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

	for (size_t i = 0; i < monomerCount(); i++)
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

double Polymer::vsRefineBackbone(void *object)
{
	Parser *parser = static_cast<Parser *>(object);
	Polymer *polymer = dynamic_cast<Polymer *>(parser);
	
	polymer->refineBackbone();
	return 0;
}

void Polymer::refineBackbone()
{
	OptionsPtr options = Options::getRuntimeOptions();
	CrystalPtr crystal = options->getActiveCrystal();
	
	const int windowSize = 10;
	const int checkSize = 5;
	RefinementType rType = RefinementFine;
	
	int skip = 1;

	while (true)
	{
		int start = (skip > 0) ? getAnchor() : getAnchor() - 1;
		int end = (skip > 0) ? monomerCount() : 0;
		std::cout << "Refining from " << start << " to " << end;
		std::cout << " using skip of " << skip << std::endl;

		for (int i = start; i != end; i += skip * windowSize)
		{
			double diff = 0;
			if (i > monomerCount() || i < 0)
			{
				break;
			}

			std::cout << "Refining using correlation with density." << std::endl;
			addParamType(ParamOptionTorsion, 0.02);
			addParamType(ParamOptionKick, 0.010);
			addParamType(ParamOptionMagicAngles, 5.0);
			addParamType(ParamOptionNumBonds, 10);
			diff += refineRange(i, i + skip * windowSize, crystal, rType);

			/*
			BoneDensity density;
			density.setCrystal(crystal);
			density.setPolymer(shared_from_this());
			density.setRange(i + skip * windowSize, i);
			density.analyse();


			BackboneState state = density.stateOfBackbone(i + skip * windowSize,
			                                              i);

			*/

			std::cout << "(Difference: " << diff << ")" << std::endl;

			if (diff < -0.008)
			{
				std::cout << "Squeezing chain to reduce expansion." << std::endl;
				refineRange(i, i + skip * windowSize,
				            crystal, RefinementRMSDZero);

				std::cout << "Re-refining torsion angles." << std::endl;
				addParamType(ParamOptionTorsion, 0.02);
				addParamType(ParamOptionNumBonds, 10);
				diff += refineRange(i, i + skip * windowSize, crystal, rType);

			}
			else if (diff != diff)
			{
				break;
			}
			else
			{
				std::cout << "Chain looking OK." << std::endl;
			}
		}

		if (skip < 0)
		{
			break;
		}
		else
		{
			skip = -1;
		}
	}

}

void Polymer::refineMonomer(MonomerPtr monomer, CrystalPtr target,
                            RefinementType rType)
{
	if (!monomer)
	{
		return;
	}

	monomer->refine(target, rType);
}

void Polymer::refineToEnd(int monNum, CrystalPtr target, RefinementType rType)
{
	int start = monNum;
	int end = (monNum < _anchorNum) ? -1 : monomerCount();
	
	refineRange(start, end, target, rType);
}

double Polymer::refineRange(int start, int end, CrystalPtr target, RefinementType rType)
{
	int skip = (start < _anchorNum) ? -1 : 1;
	if ((_anchorNum > start && _anchorNum < end) ||
	    (_anchorNum > end && _anchorNum < start))
	{
		shout_at_helen("Trying to refine a range which straddles the anchor.");	
	}

	Timer timer("refine range", true);

	int count = 0;
	double startCCAve = 0;
	std::map<MonomerPtr, double> preScores;

	std::cout << "Refining chain " << getChainID();
	std::cout  << " from residue " << start << " to ";
	std::cout << (skip > 0 ? "C" : "N");
	std::cout <<  "-terminus (residue " << end << ") ..." << std::endl;

	if (rType == RefinementModelRMSDZero)
	{
		for (int i = start; i != end; i += skip)
		{
			MonomerPtr monomer = getMonomer(i);
			if (!monomer)
			{
				continue;
			}

			for (int j = 0; j < monomer->atomCount(); j++)
			{
				AtomPtr atom = monomer->atom(j);
				atom->getModel()->getFinalPositions();	
				vec3 pos = atom->getAbsolutePosition();
				atom->setInitialPosition(pos);
			}
		}
	}

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
		if (i % 10 == 0)
		{
			std::cout << i;
		}
		
		std::cout << "\t";
		MonomerPtr monomer = getMonomer(i);
		if (!monomer)
		{
			continue;
		}

		copyParams(monomer->getSidechain());
		copyParams(monomer->getBackbone());
		refineMonomer(monomer, target, rType);
		monomer->getSidechain()->clearParams();
		monomer->getBackbone()->clearParams();

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

		std::cout << std::endl;

		endCCAve += score;
		count++;
	}

	endCCAve /= (double)count;

	std::cout << "Average CC of monomers went ";
	std::cout << ((endCCAve < startCCAve) ? "up " : "down ");
	std::cout << "from " << -startCCAve * 100 << " to " << -endCCAve * 100;
	std::cout << "." << std::endl;
	
	timer.report();

	clearParams();
	
	return -endCCAve + startCCAve;
}

void Polymer::refineVScript(void *object, RefinementType rType)
{
	OptionsPtr options = Options::getRuntimeOptions();
	CrystalPtr crystal = options->getActiveCrystal();
	
	Parser *parser = static_cast<Parser *>(object);
	Polymer *polymer = dynamic_cast<Polymer *>(parser);
	polymer->refine(crystal, rType);
	crystal->addComment("Refining against " + Options::rTypeString(rType));
}

double Polymer::vsRefinePositionsToPDB(void *object)
{
	refineVScript(object, RefinementModelPos);
	return 0;
}

double Polymer::vsRefineSidechainsToDensity(void *object)
{
	refineVScript(object, RefinementSidechain);
	return 0;	
}


void Polymer::refine(CrystalPtr target, RefinementType rType)
{
	if (rType == RefinementSidechain)
	{
		refineToEnd(getAnchor() - 1, target, rType);
		refineToEnd(getAnchor(), target, rType);
		
		return;	
	}
	
	Timer timer = Timer("refinement", true);
	
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

	timer.report();
}

std::string Polymer::makePDB(PDBType pdbType, CrystalPtr crystal)
{
	std::ostringstream stream;
	stream << getPDBContribution(pdbType, crystal);

	return stream.str();
}

void Polymer::differenceGraphs(std::string graphName, CrystalPtr diffCrystal)
{
	CSVPtr perCA = CSVPtr(new CSV(4, "resnum", "cc", "diffbackbone", "diffsidechain"));

	std::vector<double> tempCCs, tempDiffCCs, tempNs;
	double sumCC = 0;
	FFTPtr fft = diffCrystal->getFFT();
	FFTPtr difft = diffCrystal->getDiFFT();

	for (int n = 0; n < monomerCount(); n++)
	{
		if (!getMonomer(n) || !getMonomer(n)->getBackbone())
		{
			continue;
		}

		MonomerPtr monomer = getMonomer(n);

		BackbonePtr backbone = monomer->getBackbone();
		SidechainPtr sidechain = monomer->getSidechain();
		double cc = -monomer->scoreWithMap(ScoreTypeCorrel, diffCrystal);
		sumCC += cc;

		double backboneCC = -backbone->scoreWithMap(ScoreTypeMultiply,
		                                            diffCrystal, false,
		                                            MapScoreFlagDifference);
		
		double sidechainCC = -sidechain->scoreWithMap(ScoreTypeMultiply,
		                                              diffCrystal, false,
		                                              MapScoreFlagDifference);

		perCA->addEntry(4, (double)n, cc, backboneCC, sidechainCC);
	}

	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = "diffmap_" + graphName;
	plotMap["height"] = "700";
	plotMap["width"] = "1200";
	plotMap["xHeader0"] = "resnum";
	plotMap["yHeader0"] = "cc";
	plotMap["yMin0"] = "0";
	plotMap["yMax0"] = "1";

	plotMap["colour0"] = "black";
	plotMap["xTitle0"] = "Residue number";
	plotMap["yTitle0"] = "Correlation coefficient";
	plotMap["style0"] = "line";

	perCA->setSubDirectory("density_ccs");
	perCA->writeToFile("stats_" + graphName + ".csv");
	perCA->plotPNG(plotMap);
}

void Polymer::weightStrands()
{
	CSVPtr strandsCSV = CSVPtr(new CSV(1, "resnum"));
	CSVPtr sumstrandsCSV = CSVPtr(new CSV(2, "resnum", "sum"));
	
	std::vector<BondSample> positions = getAnchorModel()->getFinalPositions();

	for (int i = 0; i < positions.size(); i++)
	{
		strandsCSV->addHeader("strand_" + i_to_str(i));		
	}

	for (int j = 0; j < monomerCount(); j++)
	{
		if (!getMonomer(j))
		{
			continue;
		}

		AtomPtr ca = getMonomer(j)->findAtom("CA");
		std::vector<BondSample> positions = ca->getModel()->getFinalPositions();
		double bfac = ca->getBFactor();
		vec3 average = ca->getAbsolutePosition();
		std::vector<double> values;
		values.push_back(j);

		for (int i = 0; i < positions.size(); i++)
		{
			vec3 one = positions[i].start;
			vec3 diff = vec3_subtract_vec3(one, average);
			double length = vec3_length(diff);
			values.push_back(length);
		}
		
		strandsCSV->addEntry(values);
	}
	
	double total = 0;
	double mult = 1.5;

	for (int i = 1; i < strandsCSV->headerCount(); i++)
	{
		double sum = 0;

		for (int j = 0; j < strandsCSV->entryCount(); j++)
		{
			double value = strandsCSV->valueForEntry(i, j);
			value *= value;
			sum += value;
		}
		
		sum /= strandsCSV->entryCount();
		sum = sqrt(sum);
		
		sumstrandsCSV->addEntry(2, (double)i, sum);
		
		total += exp(-mult * sum * sum);
	}
	
	strandsCSV->writeToFile("strands.csv");
	sumstrandsCSV->writeToFile("sumstrands.csv");
	
	std::vector<double> occupancies;

	for (int i = 0; i < sumstrandsCSV->entryCount(); i++)
	{
		double sum = sumstrandsCSV->valueForEntry("sum", i);
		sum = exp(-mult * sum * sum);
		double occupancy = sum / total;
		
		printf("Occupancy %.5f\n", occupancy);

		occupancies.push_back(occupancy);
	}
	
	AbsolutePtr abs = ToAbsolutePtr(getAnchorModel());
	abs->setOccupancies(occupancies);
	propagateChange();
}

void Polymer::graph(std::string graphName)
{
	CSVPtr csv = CSVPtr(new CSV(5, "resnum", "newB", "oldB",
	                            "pos", "sidepos"));
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

			double posDisp = ca->posDisplacement();
			double sideDisp = sidechain->getAverageDisplacement();

			csv->addEntry(5, value, meanSq, ca->getInitialBFactor(),
			              posDisp, sideDisp);
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
	plotMap["filename"] = "mainchain_" + graphName;
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

	csv->setSubDirectory("bfactor_plots");
	csv->plotPNG(plotMap);
	csv->writeToFile("mainchain_" + graphName + ".csv");

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

		csv->setSubDirectory("positional_error_plots");
		sidechainCsv->setSubDirectory("bfactor_plots");
		sidechainCsv->plotPNG(plotMap);
		sidechainCsv->writeToFile("sidechain_" + graphName + ".csv");
	}

	{
		std::map<std::string, std::string> plotMap;
		plotMap["filename"] = "mainchain_" + graphName;
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

		csv->setSubDirectory("positional_error_plots");
		csv->plotPNG(plotMap);
		csv->writeToFile("mainchain_" + graphName + ".csv");
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

void Polymer::hydrogenateContents()
{
	Hydrogenator hydrogenator;

	for (int i = 0; i < monomerCount(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}	
		
		hydrogenator.setMonomer(getMonomer(i));
		hydrogenator.hydrogenate();
	}
}


double Polymer::getBackboneKick(void *object)
{
	return static_cast<Polymer *>(object)->_kick;
}

void Polymer::setBackboneKick(void *object, double value)
{
	Polymer *poly = static_cast<Polymer *>(object);
	int anchor = poly->getAnchor();
	poly->_kick = value;

	for (int i = 0; i < poly->atomCount(); i++)
	{
		ModelPtr model = poly->atom(i)->getModel();
		if (!model || !model->isBond())
		{
			continue;
		}
		
		BondPtr bond = ToBondPtr(model);
		
		if (!bond->connectsAtom("CA"))
		{
			continue;	
		}
		
		double mult = 1;
		if (poly->atom(i)->getResidueNum() < anchor)
		{
			mult = -1;
		}
		
		Bond::setTorsionBlur(&*bond, value * mult);
	}
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

	for (size_t i = 0; i < finals->size(); i++)
	{
		sum = vec3_add_vec3(sum, finals->at(i).start);
	}

	vec3_mult(&sum, 1/(double)finals->size());
	vec3 total = empty_vec3();

	for (size_t i = 0; i < finals->size(); i++)
	{
		vec3 onePos = finals->at(i).start;
		
		if (_centroidOffsets.size() > i)
		{
			vec3 offset = _centroidOffsets[i];
			vec3_mult(&offset, -1);
			vec3_add_to_vec3(&onePos, offset);
		}
		
		vec3_add_to_vec3(&total, onePos);
		
		results.push_back(onePos);
	}
	
	vec3_mult(&total, 1 / (double)finals->size());
	
	for (size_t i = 0; i < results.size(); i++)
	{
		results[i] = vec3_subtract_vec3(results[i], total);
	}

	return results;
}

void Polymer::applyPolymerChanges()
{
	for (size_t i = 0; i < atomCount(); i++)
	{
		if (!atom(i)) continue;

		ModelPtr model = atom(i)->getModel();
		model->setPolymerChanged();	
	}
}

void Polymer::applyTranslationTensor()
{
	_transTensorOffsets.clear();
	_extraRotationMats.clear();

	std::vector<vec3> sphereDiffs = getAnchorSphereDiffs();

	vec3 sum = empty_vec3();
	double nonExpLength = 0;
	double expLength = 0;

	for (size_t i = 0; i < sphereDiffs.size(); i++)
	{
		vec3 diff = sphereDiffs[i];
		vec3_add_to_vec3(&diff, _sphereDiffOffset);
		vec3 diffTensored = diff;
		mat3x3_mult_vec(_transTensor, &diffTensored);
		vec3 movement = vec3_subtract_vec3(diffTensored, diff);
		
		double length = vec3_length(movement);

		double mult = exp(length * _transExponent);
		vec3_mult(&movement, _overallScale);

		nonExpLength += vec3_length(movement);

		vec3_mult(&movement, mult);
		expLength += vec3_length(movement);

		_transTensorOffsets.push_back(movement);
		vec3_add_to_vec3(&sum, movement);
	}
	
	double normalise = nonExpLength / expLength;
	
	if (normalise != normalise)
	{
		normalise = 1;
	}
	
	vec3_mult(&sum, -1 / (double)_transTensorOffsets.size());

	for (size_t i = 0; i < _transTensorOffsets.size(); i++)
	{
		vec3_mult(&_transTensorOffsets[i], normalise);
//		vec3_add_to_vec3(&_transTensorOffsets[i], sum);
	}

	applyPolymerChanges();
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

	for (size_t i = 0; i < sphereDiffs.size(); i++)
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
	
	applyPolymerChanges();
}

void Polymer::superimpose()
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	crystal->addComment("Superimposing ensemble for chain " + getChainID());
	
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

		for (size_t i = 0; i < poss.size(); i++)
		{
			vec3 pos = poss[i];
			csv->addEntry(3, pos.x, pos.y, pos.z);
		}

//		csv->writeToFile("ca_pos.csv");
	}

	if (model->isAbsolute())
	{
		std::vector<vec3> sphereAngles = ToAbsolutePtr(model)->getSphereAngles();
		CSVPtr csv = CSVPtr(new CSV(10, "psi", "phi", "theta", "corr_x",
		                            "corr_y", "corr_z", "rot_angle",
		                            "rot_axis_x", "rot_axis_y", "rot_axis_z"));
		CSVPtr one = CSVPtr(new CSV(3, "x", "y", "z"));

		for (size_t i = 0; i < sphereAngles.size(); i++)
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

//		csv->writeToFile("kabsch.csv");
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

	double totalWeights = 0;
	getAnchorModel()->getFinalPositions();
	vec3 oldPos = getAnchorModel()->getAbsolutePosition();

	/* For every cAlpha atom in the structure, get the mean centroid
	 * position of each conformer... */
	for (size_t i = 0; i < monomerCount(); i++)
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

		/* Prepare the vectors for the first time if necessary */
		if (!addedVecs.size())
		{
			addedVecs = std::vector<vec3>(samples->size(), make_vec3(0, 0, 0));;
		}

		double weight = 1;//ca->getBFactor();
		totalWeights += weight;

		/* Sum over the uncorrected start positions */
		for (size_t j = 0; j < samples->size(); j++)
		{
			vec3 added = samples->at(j).start;

			vec3_mult(&added, weight);

			addedVecs[j] = vec3_add_vec3(addedVecs[j], added);
		}
	}

	/* Divide each by total weights */
	double inverse = 1 / totalWeights;

	for (size_t j = 0; j < addedVecs.size(); j++)
	{
		vec3_mult(&addedVecs[j], inverse);
	}

	/* Mean of all the centroids in the structure */
	vec3 meanPos = make_vec3(0, 0, 0);

	for (size_t i = 0; i < addedVecs.size(); i++)
	{
		vec3 addition = addedVecs[i];
		meanPos = vec3_add_vec3(meanPos, addition);
	}

	vec3_mult(&meanPos, 1 / (double)addedVecs.size());

	// starting point for rotation centre of whole molecule movements
	// Used by Model to correct for translation stuff
	_rotationCentre = meanPos;

	// Remove the old centroid positions if not already
	_centroidOffsets.clear();

	// Find the offsets to bring all centroids to the mean value
	for (size_t i = 0; i < addedVecs.size(); i++)
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

	for (size_t i = 0; i < _centroidOffsets.size(); i++)
	{
		vec3 offset = vec3_add_vec3(_centroidOffsets[i], backToOld);
		_centroidOffsets[i] = offset;
	}

	std::cout << "Calculated corrections for " << _centroidOffsets.size() << " structures." << std::endl;
}

void Polymer::minimiseRotations()
{
	size_t num = 0;

	std::vector<vec3> tmpCentroids;

	AtomGroupPtr backbone = getAllBackbone();

	/* Find the number of samples */
	for (size_t i = 0; i < backbone->atomCount(); i++)
	{
		AtomPtr atom = backbone->atom(i);

		std::vector<BondSample> *samples;
		samples = atom->getModel()->getManyPositions();

		num = samples->size();
		break;
	}
	
	std::vector<mat3x3> tmpMats;

	for (size_t i = 0; i < num; i++)
	{
		std::vector<double> weights;
		std::vector<vec3> fixedVecs, variantVecs;

		for (size_t k = 0; k < backbone->atomCount(); k++)
		{
			AtomPtr anAtom = backbone->atom(k);
			if (!anAtom) continue;

			std::vector<BondSample> *samples;
			ModelPtr model = anAtom->getModel();

			if (!model)
			{
				shout_at_helen("Missing model for backbone atom!");
			}

			double weight = 1;//anAtom->getBFactor();
			weights.push_back(weight);

			samples = model->getManyPositions();
			vec3 fixed = model->getAbsolutePosition();

			vec3 variant = samples->at(i).start;
			fixedVecs.push_back(fixed);
			variantVecs.push_back(variant);
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

	for (size_t i = 0; i < atomCount(); i++)
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

	for (size_t i = 0; i < atomCount(); i++)
	{
		if (!atom(i)->getModel())
		{
			std::cout << "Missing model for " << atom(i)->shortDesc() <<
			std::endl;
			bondsOk = 0;
		}

		if (atom(i)->getModel()->isBond())
		{
			bondsOk *= ToBondPtr(atom(i)->getModel())->test();
		}
	}

	if (!bondsOk)
	{
		std::cout << "Bonds FAILED test for polymer " << getChainID() <<
		std::endl;
	}

	return true;
}

void Polymer::reportParameters()
{
	int count = 0;

	for (size_t i = 0; i < atomCount(); i++)
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

void Polymer::vsTransTensorOverall(void *object, double value)
{
	Parser *parser = static_cast<Parser *>(object);
	Polymer *polymer = dynamic_cast<Polymer *>(parser);

	polymer->_transTensor = make_mat3x3();
	mat3x3_mult_scalar(&polymer->_transTensor, value);
	polymer->applyTranslationTensor();

	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	crystal->addComment("Changing overall translation scale to "
	                    + f_to_str(value, 2) + " for chain " + 
	                    polymer->getChainID());
}

double Polymer::vsFitTranslation(void *object)
{
	Parser *parser = static_cast<Parser *>(object);
	Polymer *polymer = dynamic_cast<Polymer *>(parser);

	polymer->optimiseWholeMolecule(true, false);
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	crystal->addComment("Refining overall translation matrix for chain "
	                    + polymer->getChainID());
}

double Polymer::vsFitRotation(void *object)
{
	Parser *parser = static_cast<Parser *>(object);
	Polymer *polymer = dynamic_cast<Polymer *>(parser);

	polymer->optimiseWholeMolecule(false, true);
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	crystal->addComment("Refining overall rotation parameters for chain "
	                    + polymer->getChainID());
}

void Polymer::optimiseWholeMolecule(bool translation, bool rotation)
{
	std::cout << "Optimising whole molecule shifts to match the electron density." << std::endl;

	Timer timer("whole molecule fit", true);
	
	FlexGlobal target;
	NelderMeadPtr nelderMead = NelderMeadPtr(new NelderMead());

	if (translation)
	{
		attachTargetToRefinement(nelderMead, target);
		nelderMead->addParameter(this, getOverallScale, setOverallScale,
		                         0.5, 0.01, "overall_scale");
		
//		nelderMead->addParameter(this, getTransExponent, setTransExponent, 0.5, 0.01);
		nelderMead->setCycles(25);
		nelderMead->setVerbose(true);
		nelderMead->refine();
		nelderMead->clearParameters();
	}
	
	attachTargetToRefinement(nelderMead, target);

	if (translation)
	{
		nelderMead->addParameter(this, getTransTensor11, setTransTensor11,
		                         0.5, 0.01, "t11");
		nelderMead->addParameter(this, getTransTensor12, setTransTensor12,
		                         0.1, 0.01, "t12");
		nelderMead->addParameter(this, getTransTensor21, setTransTensor21,
		                         0.1, 0.01, "t21");
		nelderMead->addParameter(this, getTransTensor13, setTransTensor13,
		                         0.1, 0.01, "t13");
		nelderMead->addParameter(this, getTransTensor22, setTransTensor22,
		                         0.5, 0.01, "t22");
		nelderMead->addParameter(this, getTransTensor31, setTransTensor31,
		                         0.1, 0.01, "t31");
		nelderMead->addParameter(this, getTransTensor23, setTransTensor23,
		                         0.1, 0.01, "t23");
		nelderMead->addParameter(this, getTransTensor32, setTransTensor32,
		                         0.1, 0.01, "t32");
		nelderMead->addParameter(this, getTransTensor33, setTransTensor33,
		                         0.5, 0.01, "t33");
	}

	if (rotation)
	{
		nelderMead->addParameter(this, getRotPhi, setRotPhi, 0.2, 0.0001);
		nelderMead->addParameter(this, getRotPsi, setRotPsi, 0.2, 0.0001);
		nelderMead->addParameter(this, getRotAngle, setRotAngle, 0.002, 0.0001);
		nelderMead->addParameter(this, getRotExponent, setRotExponent, 0.5, 0.01);
		nelderMead->addParameter(this, getRotCentreX, setRotCentreX, 2.0, 0.01);
		nelderMead->addParameter(this, getRotCentreY, setRotCentreY, 2.0, 0.01);
		nelderMead->addParameter(this, getRotCentreZ, setRotCentreZ, 2.0, 0.01);
	}

	nelderMead->setCycles(25);
	nelderMead->setVerbose(true);

	nelderMead->refine();
	nelderMead->clearParameters();

	timer.report();
}

AtomGroupPtr Polymer::getAllBackbone()
{
	if (_allBackbones)
	{
		return _allBackbones;
	}
	
	_allBackbones = AtomGroupPtr(new AtomGroup());

	for (size_t i = 0; i < monomerCount(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}
		
		BackbonePtr bone = getMonomer(i)->getBackbone();
		if (bone)
		{
			_allBackbones->addAtomsFrom(bone);
		}
	}
	
	return _allBackbones;
}

void Polymer::attachTargetToRefinement(RefinementStrategyPtr strategy,
                                       FlexGlobal &target)
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	AtomGroupPtr allBackbone = getAllBackbone();
	target.setAtomGroup(allBackbone);
	target.setCrystal(crystal);
	target.matchElectronDensity();
	strategy->setEvaluationFunction(FlexGlobal::score, &target);
	FlexGlobal::score(&target);
}

double Polymer::vsFindKickAndDampen(void *object)
{
	Parser *parser = static_cast<Parser *>(object);
	Polymer *polymer = dynamic_cast<Polymer *>(parser);
	
	return findOverallKickAndDampen(polymer);
}

double Polymer::findOverallKickAndDampen(void *object)
{
	Polymer *poly = static_cast<Polymer *>(object);
	NelderMeadPtr nelderMead = NelderMeadPtr(new NelderMead());

	nelderMead->addParameter(poly, getBackboneKick, setBackboneKick,
	                         0.001, 0.0001);
	nelderMead->addParameter(poly, getBackboneDampening, setBackboneDampening,
	                         0.005, 0.0005);
	nelderMead->setVerbose(true);
	
	Timer timer("overall kick and dampen", true);

	FlexGlobal target;
	poly->attachTargetToRefinement(nelderMead, target);
	nelderMead->refine();

	timer.report();

	return 0;
}

double Polymer::vsSandbox(void *object)
{
	Parser *parser = static_cast<Parser *>(object);
	Polymer *polymer = dynamic_cast<Polymer *>(parser);
	
	polymer->refineLoop(67, false);
	polymer->refineLoop(66, false);
	polymer->refineLoop(67, true);
	polymer->refineLoop(66, true);
	return 0;
	
	int anchorNum = polymer->getAnchor();
	polymer->refineEverything(anchorNum - 2);
	polymer->refineEverything(anchorNum + 2);
}

void Polymer::refineLoop(int start, bool magic)
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	MonomerPtr monomer = getMonomer(start);
	setupNelderMead();
	setCrystal(crystal);

	if (!magic)
	{
		addParamType(ParamOptionKick, 0.005);
	}
	else
	{
		addParamType(ParamOptionMagicAngles, 20);
	}
	
	if (!monomer) return;
	AtomPtr nitro = monomer->findAtom("CA");
	if (!nitro) return;
	BondPtr nitroBond = ToBondPtr(nitro->getModel());
	if (!nitroBond) return;
	setupTorsionSet(nitroBond, 0, 300, 0, 0, 0, 0);

	setScoreType(ScoreTypeBFactorAgreement);
	
//	setSilent();
	setCycles(200);
	
	sample();
	clearParams();
}

void Polymer::refineEverything(int start)
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	MonomerPtr monomer = getMonomer(start);
	
	if (!monomer) return;
	
	AtomPtr nitro = monomer->findAtom("N");
	
	if (!nitro) return;
	
	BondPtr nitroBond = ToBondPtr(nitro->getModel());
	
	if (!nitroBond) return;
	
	setupNelderMead();
	setCrystal(crystal);
	
	addParamType(ParamOptionMagicAngles, 20);
	/*
	addParamType(ParamOptionKick, 0.05);
	*/
	setupTorsionSet(nitroBond, 0, 1000, 0, 0, 0, 0);
	
	setVerbose();
	setCycles(100);
	
	sample();
	clearParams();
	return;

	setupNelderMead();
	setCrystal(crystal);
	
	addParamType(ParamOptionMagicAngles, 20);
	setupTorsionSet(nitroBond, 0, 1000, 0, 0, 0, 0);
	
	setVerbose();
	setCycles(400);
	
	sample();
	clearParams();
}

void Polymer::addProperties()
{
	Molecule::addProperties();

	addIntProperty("anchor_res", &_anchorNum);
	addMat3x3Property("trans_tensor", &_transTensor);
	addDoubleProperty("overall_scale", &_overallScale);

	for (size_t i = 0; i < monomerCount(); i++)
	{
		if (!getMonomer(i)) continue;

		addChild("monomer", getMonomer(i));
	}
	
	exposeFunction("set_rot_angle", vsSetRotAngle);
	exposeFunction("refine_positions_to_pdb", vsRefinePositionsToPDB);
	exposeFunction("refine_sidechains_to_density", vsRefineSidechainsToDensity);
	exposeFunction("refine_backbone_to_density", vsRefineBackbone);
	exposeFunction("superimpose", vsSuperimpose);
	exposeFunction("set_overall_translation", vsTransTensorOverall);
	exposeFunction("fit_translation", vsFitTranslation);
	exposeFunction("fit_rotation", vsFitRotation);
	exposeFunction("sandbox", vsSandbox);
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
