//
//  Polymer.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Polymer.h"
#include "Timer.h"
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
#include "Bond.h"
#include "CSV.h"
#include <sstream>
#include <iomanip>
#include "maths.h"
#include <float.h>
#include "FileReader.h"
#include "Options.h"
#include "FlexGlobal.h"
#include "FlexLocal.h"
#include "Reflex.h"
#include "RefinementNelderMead.h"
#include "RefinementGridSearch.h"
#include "Hydrogenator.h"

Polymer::Polymer()
{
	_kick = Options::getKick();
	_kickShift = 0.1;
	_anchorNum = -1;
	_totalMonomers = 0;
	_startB = Options::getBStart();
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

void Polymer::reflex()
{
	Reflex reflex;
	reflex.setPolymer(shared_from_this());
	reflex.setPieceCount(1);
	reflex.calculate();
}

void Polymer::refineLocalFlexibility()
{
	FlexLocal local;
	local.setPolymer(shared_from_this(), _kickShift);
	local.refine();
	_kickShift = local.getShift();
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
	AtomPtr to_c = getMonomer(_anchorNum)->findAtom("CA");
	AtomPtr to_n = getMonomer(_anchorNum - 1)->findAtom("C");

	/* Specify heavy alignment atoms around the anchor point */
	AtomPtr ca = getMonomer(_anchorNum)->findAtom("CA");
	AtomPtr c = getMonomer(_anchorNum)->findAtom("C");
	AtomPtr prev_c = getMonomer(_anchorNum - 1)->findAtom("C");
	AtomPtr prev_ca = getMonomer(_anchorNum - 1)->findAtom("CA");

	ModelPtr nModel = n->getModel();

	if (nModel->isAbsolute())
	{
		AnchorPtr newAnchor = AnchorPtr(new Anchor(ToAbsolutePtr(nModel)));
		newAnchor->setBFactor(_startB);
		newAnchor->setNeighbouringAtoms(prev_ca, to_n, to_c, c);
		n->setModel(newAnchor);
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
	
	BondPtr n2ca = ToBondPtr(ca->getModel());
	BondPtr ca2c = ToBondPtr(c->getModel());
	BondPtr n2c = ToBondPtr(prev_c->getModel());
	BondPtr c2ca = ToBondPtr(prev_ca->getModel());
	
	ca2c->setHeavyAlign(prev_c);
	n2ca->setHeavyAlign(prev_ca);
	n2c->setHeavyAlign(c);
	c2ca->setHeavyAlign(ca);
	
	n2ca->checkForSplits(shared_from_this());
	n2c->checkForSplits(shared_from_this());
}

void Polymer::removeAtom(AtomPtr atom)
{
	MonomerPtr mon = atom->getMonomer();
	mon->removeAtom(atom);
	
	AtomGroup::removeAtom(atom);
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

	std::cout << "| I am a polymer with " << _monomers.size() << " monomers."
	<< std::endl;
}

double Polymer::vsRefineBackbone(void *object)
{
	Parser *parser = static_cast<Parser *>(object);
	Polymer *polymer = dynamic_cast<Polymer *>(parser);
	
	polymer->refineBackbone();
	return 0;
}

void Polymer::vsRefineBackboneFrom(void *object, double position)
{
	Parser *parser = static_cast<Parser *>(object);
	Polymer *polymer = dynamic_cast<Polymer *>(parser);
	
	polymer->refineBackboneFrom(position);
}

void Polymer::refineBackboneFrom(int position)
{
	OptionsPtr options = Options::getRuntimeOptions();
	CrystalPtr crystal = options->getActiveCrystal();
	
	RefinementType rType = RefinementFine;
	
	std::cout << "Refining magic angles using correlation with density."
	 << std::endl;

	addParamType(ParamOptionMagicAngles, 36.0);
	addParamType(ParamOptionNumBonds, 4);

	refineToEnd(position, crystal, rType);
}

void Polymer::refineBackbone()
{
	int anchor = getAnchor();

	refineBackboneFrom(anchor - 5);
	refineBackboneFrom(anchor + 5);
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
	std::map<MonomerPtr, vec3> preScores;

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
				atom->getModel()->refreshPositions();	
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
		
		BackbonePtr bone = monomer->getBackbone();
		double backScore = bone->scoreWithMap(ScoreTypeCorrel, target);
		SidechainPtr side = monomer->getSidechain();
		double sideScore = side->scoreWithMap(ScoreTypeCorrel, target);

		vec3 scores = make_vec3(score, backScore, sideScore);
		preScores[monomer] = scores;
		count++;
	}

	startCCAve /= (double)count;
	count = 0;

	double endCCAve = 0;
	
	std::cout << "\t  Backbone             | Sidechain" << std::endl;
	bool changed = true;

	for (int i = start; i != end; i += skip)
	{
		if (i % 10 == 0 && changed == true)
		{
			std::cout << i;
			changed = false;
		}

		std::cout << "\t";
		MonomerPtr monomer = getMonomer(i);
		
		if (!monomer) continue;
		
		BackbonePtr bone = monomer->getBackbone();
		SidechainPtr side = monomer->getSidechain();

		if (!monomer)
		{
			continue;
		}

		changed = true;
		
		copyParams(side);
		copyParams(bone);
		refineMonomer(monomer, target, rType);
		side->clearParams();
		bone->clearParams();

		double score = monomer->scoreWithMap(ScoreTypeCorrel, target);	
		double backScore = bone->scoreWithMap(ScoreTypeCorrel, target);	
		double sideScore = side->scoreWithMap(ScoreTypeCorrel, target);	
		double pre = preScores[monomer].x;
		double preBack = preScores[monomer].y;
		double preSide = preScores[monomer].z;

		/* Scores are negative by default */
		double diff = (score - pre) * 100.;
		double backdiff = (backScore - preBack) * 100.;
		double sidediff = (sideScore - preSide) * 100.;

		/* improvement of 2% will be 20 + symbols after the residue! */
		/* display backbone on one side */
		int signs = fabs(backdiff * 10);
		int dir = (backdiff < 0);	
		
		if (signs > 20) signs = 20;

		std::cout << " ";
		for (int j = 0; j < signs; j++)
		{
			std::cout << (dir ? "+" : "-");
		}
		for (int j = signs; j < 20; j++)
		{
			std::cout << " ";	
		}

		signs = fabs(sidediff * 10);
		dir = (sidediff < 0);	

		std::cout << " | ";
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
	target->addComment("Refining against " + Options::rTypeString(rType));

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

void Polymer::graph(std::string graphName)
{
	ramachandranPlot();

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
			caBlur = Bond::getKick(&*caBond);
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
	
	bool extra = Options::makeDiagnostics();

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

	if (extra)
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

	if (extra)
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

void Polymer::ramachandranPlot()
{
	CSVPtr csv = CSVPtr(new CSV(3, "res", "phi", "psi"));

	for (int i = 0; i < monomerCount(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}
		
		bool forwards = (i > getAnchor());

		AtomList phiAtoms = getMonomer(i)->findAtoms("N");
		AtomList psiAtoms = getMonomer(i)->findAtoms("C");
		
		if (phiAtoms.size() != psiAtoms.size())
		{
			continue;
		}

		for (int j = 0; j < phiAtoms.size(); j++)
		{
			ModelPtr mPhi = phiAtoms[j]->getModel();
			ModelPtr mPsi = psiAtoms[j]->getModel();
			
			if (!mPhi->isBond() || !mPsi->isBond())
			{
				continue;
			}
			
			double tPhi = Bond::getTorsion(&*ToBondPtr(mPhi));
			double tPsi = Bond::getTorsion(&*ToBondPtr(mPsi));
			
			csv->addEntry(3, i, rad2deg(tPhi), rad2deg(tPsi));
		}
	}

	csv->writeToFile("ramachandran.csv");
}

std::string Polymer::makePDB(PDBType pdbType, CrystalPtr crystal, 
                             int conformer)
{
	std::ostringstream stream;
	for (int i = 0; i < monomerCount(); i++)
	{
		MonomerPtr monomer = getMonomer(i);
		
		if (!monomer)
		{
			continue;
		}

		stream << monomer->getPDBContribution(pdbType, crystal, conformer);
	}

	return stream.str();
	
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

		n->getModel()->refreshPositions();
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

void Polymer::vsMultiplyBackboneKick(void *object, double value)
{
	Parser *parser = static_cast<Parser *>(object);
	Polymer *poly = dynamic_cast<Polymer *>(parser);

	for (int i = 0; i < poly->monomerCount(); i++)
	{
		MonomerPtr mono = poly->getMonomer(i);
		
		if (!mono) continue;

		BackbonePtr back = mono->getBackbone();
		
		for (int i = 0; i < back->atomCount(); i++)
		{
			AtomPtr atom = back->atom(i);
			ModelPtr model = atom->getModel();
			if (!model || !model->isBond())
			{
				continue;
			}

			BondPtr bond = ToBondPtr(model);
			double kick = Bond::getKick(&*bond);
			kick *= value;
			Bond::setKick(&*bond, kick);
		}
	}
	
	poly->propagateChange();
	poly->refreshPositions();
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
				Bond::setKick(&*bond, kick);
			}
		}
	}
}

double Polymer::getInitialKick(void *object)
{
	Polymer *polymer = static_cast<Polymer *>(object);
	int monomerNum = polymer->_anchorNum;
	return polymer->getMonomer(monomerNum)->getKick();
}

ExplicitModelPtr Polymer::getAnchorModel()
{
	MonomerPtr anchoredRes = getMonomer(getAnchor());
	if (!anchoredRes)
	{
		return ExplicitModelPtr();
	}

	ModelPtr model = anchoredRes->findAtom("N")->getModel();

	return ToExplicitModelPtr(model);
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

double Polymer::vsRefineLocalFlexibility(void *object)
{
	Parser *parser = static_cast<Parser *>(object);
	Polymer *polymer = dynamic_cast<Polymer *>(parser);
	polymer->refineLocalFlexibility();
}

void Polymer::vsOmitResidues(void *object, double start, double end)
{
	Parser *parser = static_cast<Parser *>(object);
	Polymer *polymer = dynamic_cast<Polymer *>(parser);

	polymer->downWeightResidues(start, end, 0);
}

void Polymer::vsUnomitResidues(void *object, double start, double end)
{
	Parser *parser = static_cast<Parser *>(object);
	Polymer *polymer = dynamic_cast<Polymer *>(parser);

	polymer->downWeightResidues(start, end, 1);
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
	<< start << "-" << end << " to weighting of " << value << std::endl;
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
	int torsionCountBackbone = 0;
	int angleCountBackbone = 0;
	
	int torsionCountSidechain = 0;
	int angleCountSidechain = 0;

	for (size_t i = 0; i < atomCount(); i++)
	{
		AtomPtr a = atom(i);
		bool back = false;

		if (a->isBackbone() || a->isBackboneAndSidechain())
		{
			back = true;
		}

		if (a->getModel()->isBond())
		{
			if (ToBondPtr(a->getModel())->isTorsionRefinable())
			{
				if (back)
				{
					torsionCountBackbone++;
				}
				else
				{
					torsionCountSidechain++;
				}
			}

			if (ToBondPtr(a->getModel())->getRefineBondAngle())
			{
				if (back)
				{
					angleCountBackbone++;
				}
				else
				{
					angleCountSidechain++;
				}
			}
		}
	}
	
	std::cout << std::endl;
	std::cout << "|----------------" << std::endl;
	std::cout << "| Bond-based parameter groups (" << getChainID() << "):" << std::endl;
	std::cout << "|----------------" << std::endl;
	std::cout << "|          |   Backbone |  Sidechain |" << std::endl;
	std::cout << "| Torsion  |  " << std::setw(9) << torsionCountBackbone;
	std::cout << " |  " << std::setw(9) << torsionCountSidechain;
	std::cout << " |" << std::endl;
	std::cout << "| Angle    |  " << std::setw(9) << angleCountBackbone;
	std::cout << " |  " << std::setw(9) << angleCountSidechain;
	std::cout << " |" << std::endl;
	std::cout << "|----------------" << std::endl;
	std::cout << std::endl;
}

void Polymer::refineAnchorMovements()
{
	std::cout << "Optimising anchor shifts to match the electron density." << std::endl;

	Timer timer("anchor fit", true);
	
	AnchorPtr anchor = ToAnchorPtr(getAnchorModel());
	
	FlexGlobal target;
	NelderMeadPtr nelderMead = NelderMeadPtr(new RefinementNelderMead());

	attachTargetToRefinement(nelderMead, target);

	nelderMead->setCycles(24);

	nelderMead->clearParameters();
	
	anchor->addTranslationParameters(nelderMead);
	nelderMead->refine();
	nelderMead->clearParameters();

	anchor->addLibrationParameters(nelderMead);
	nelderMead->refine();
	return;

	FlexLocal local;
	local.setPolymer(shared_from_this(), _kickShift);
	local.refineAnchor();

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
                                       FlexGlobal &target, bool isotropy)
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	AtomGroupPtr allBackbone = getAllBackbone();
	target.setAtomGroup(allBackbone);
	target.setCrystal(crystal);
	
	if (!isotropy)
	{
		target.matchElectronDensity();
	}
	else
	{
		target.maximiseIsotropy();
	}

	strategy->setEvaluationFunction(FlexGlobal::score, &target);
	FlexGlobal::score(&target);
}

double Polymer::vsSandbox(void *object)
{
	Parser *parser = static_cast<Parser *>(object);
	Polymer *polymer = dynamic_cast<Polymer *>(parser);
	
	int anchorNum = polymer->getAnchor();
}

void Polymer::addProperties()
{
	Molecule::addProperties();

	addIntProperty("anchor_res", &_anchorNum);

	for (size_t i = 0; i < monomerCount(); i++)
	{
		if (!getMonomer(i)) continue;

		addChild("monomer", getMonomer(i));
	}
	
	exposeFunction("refine_positions_to_pdb", vsRefinePositionsToPDB);
	exposeFunction("refine_sidechains_to_density", vsRefineSidechainsToDensity);
	exposeFunction("refine_local_flexibility", vsRefineLocalFlexibility);
	exposeFunction("refine_backbone_to_density", vsRefineBackbone);
	exposeFunction("refine_backbone_to_density_from", vsRefineBackboneFrom);
	exposeFunction("sandbox", vsSandbox);

	exposeFunction("omit_residues", vsOmitResidues);
	exposeFunction("unomit_residues", vsUnomitResidues);
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
	//applyTranslationTensor();
	
}


