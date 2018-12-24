//
//  Polymer.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Polymer.h"
#include "Timer.h"
#include "Twist.h"
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
#include "Refitter.h"
#include "RefinementNelderMead.h"
#include "RefinementGridSearch.h"
#include "Hydrogenator.h"

Polymer::Polymer()
{
	_kick = Options::getKick();
	_kickShift = 0.05;
	_anchorNum = -1;
	_totalMonomers = 0;
	_startB = Options::getBStart();
	_flexibilityParams = 0;
	_positionalParams = 0;
	_whacked = -1;
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
	reflex.setPieceCount(3);
	reflex.calculate();
}

void Polymer::refineLocalFlexibility()
{
	Timer timer("flexibility refinement", true);
	FlexLocal local;
	whack();
	local.setPolymer(shared_from_this(), _kickShift);
	
	if (getAnchorModel()->whackCount() > 0)
	{
		local.setWhacking(true);
	}
	
	local.refine();
	_kickShift = local.getShift();
	timer.report();
}

bool Polymer::isWhacking()
{
	return (getAnchorModel()->whackCount() > 0);
}

void Polymer::whackMonomer(MonomerPtr mon)
{
	if (!mon)
	{
		return;
	}
	
	AtomList atoms = mon->getBackbone()->topLevelAtoms();
	atoms = mon->findAtoms("CA");
	AnchorPtr anchor = getAnchorModel();
	long num = mon->getResidueNum();

	if (atoms.size() == 0)
	{
		return;
	}

	AtomPtr atom = atoms[0];

	if (!atom->getModel()->isBond())
	{
		return;
	}

	BondPtr bond = ToBondPtr(atom->getModel());

	if (bond->getMajor()->getModel()->isAnchor())
	{
		return;
	}

	WhackPtr whack = WhackPtr(new Whack());
	whack->setBond(bond);

	if (!whack->isValid())
	{
		return;
	}

	whack->addToAnchor(anchor);

	BondPtr next;

	TwistPtr twist = TwistPtr(new Twist());
	twist->setBond(bond);
	twist->addToAnchor(anchor);

	while (true)
	{
		if (bond->downstreamBondGroupCount() && 
		    bond->downstreamBondCount(0))
		{
			next = bond->downstreamBond(0, 0);
		}
		else
		{
			break;
		}

		TwistPtr twist = TwistPtr(new Twist());
		twist->setBond(next);
		twist->addToAnchor(anchor);

		if (next->getAtom()->getResidueNum() != num)
		{
			break;
		}

		next->setRefineFlexibility(false);
		bond = ToBondPtr(next);
	}

}

void Polymer::whack()
{
	if (isWhacking())
	{
		return;
	}

	std::cout << "Whacking chain " << getChainID() << std::endl;
	
	for (int i = getAnchor(); i < monomerEnd(); i++)
	{
		whackMonomer(getMonomer(i));
	}

	for (int i = getAnchor() - 1; i >= monomerBegin(); i--)
	{
		whackMonomer(getMonomer(i));
	}

	getAnchorModel()->propagateChange(-1, true);
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
	
	whack();
}

void Polymer::removeAtom(AtomPtr atom)
{
	MonomerPtr mon = atom->getMonomer();
	
	if (mon)
	{
		mon->removeAtom(atom);
	}
	
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

void refineLeftRegion(AtomGroupPtr region, CrystalPtr target, double light)
{
	region->addParamType(ParamOptionTorsion, light);
	region->addParamType(ParamOptionBondAngle, light);
	region->addParamType(ParamOptionNumBonds, 4);
	region->refine(target, RefinementSavedPos);
}

void Polymer::refineAroundMonomer(int central, CrystalPtr target)
{
	int pad = 2;
	MonomerPtr monomer = getMonomer(central);
	
	if (!monomer)
	{
		return;
	}
	
	int coreStart = central - pad;
	int coreEnd = central + pad;
	
	AtomGroupPtr coreRegion = monomerRange(coreStart, coreEnd);
	
	if (!coreRegion || coreRegion->atomCount() == 0)
	{
		return;
	}
	
	int anchor = getAnchor();
	
	double step = 1.5;

	std::cout << "Refining core region, " << getChainID() << " " << coreStart
	<< " to " << coreEnd << std::flush;

	coreRegion->addParamType(ParamOptionNumBonds, (2 * pad + 1) * 3);
	coreRegion->addParamType(ParamOptionTorsion, step);
	coreRegion->addParamType(ParamOptionTwist, step);
	coreRegion->addParamType(ParamOptionBondAngle, step);
	coreRegion->addParamType(ParamOptionMaxTries, 1);
	coreRegion->addIncludeForRefinement(coreRegion);
	coreRegion->refine(target, RefinementCrude);

	coreRegion->saveAtomPositions();
	clearTwists();
	
	bool coversAnchor = (anchor >= coreStart && anchor <= coreEnd);
	
	if (coversAnchor)
	{
		refineAnchorPosition(target);
	}
	
	std::cout << "." << std::flush;
	
	AtomGroupPtr leftRegion = monomerRange(-1, anchor - 1);
	AtomGroupPtr rightRegion = monomerRange(anchor + 1, -1);

	refineLeftRegion(leftRegion, target, step);
	refineLeftRegion(rightRegion, target, step);

	getAnchorModel()->propagateChange(-1, true);

	std::cout << " Displacement now " << getAverageDisplacement()
	<< " Å." << std::endl;
}

void Polymer::refineToEnd(int monNum, CrystalPtr target, RefinementType rType)
{
	int start = monNum;
	int end = (monNum < _anchorNum) ? -1 : monomerCount();
	
	refineRange(start, end, target, rType);
}

AtomGroupPtr Polymer::monomerRange(int start, int end)
{
	AtomGroupPtr all = AtomGroupPtr(new AtomGroup());
	
	if (start == -1)
	{
		start = monomerBegin();
	}

	if (end == -1)
	{
		end = monomerEnd();
	}
	
	for (int i = start; i <= end; i++)
	{
		MonomerPtr monomer = getMonomer(i);

		if (!monomer)
		{
			continue;
		}
		
		all->addAtomsFrom(monomer);
	}
	
	return all;
}

double Polymer::refineRange(int start, int end, CrystalPtr target, 
                            RefinementType rType)
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

		if (rType == RefinementCrude)
		{
			vec3 centre = monomer->centroid();
			Options::getRuntimeOptions()->focusOnPosition(centre);
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

void Polymer::clearTwists()
{
	getAnchorModel()->clearTwists();
}

void Polymer::refineAnchorPosition(CrystalPtr target)
{
	bool changed = true;
	int count = 0;
	
	while (changed && count < 50)
	{
		count++;
		setupNelderMead();
		setCrystal(target);
		setCycles(16);
		setScoreType(ScoreTypeSavedPos);
		setSilent(true);
		addAnchorParams(getAnchorModel());

		changed = sample();
	}

	getAnchorModel()->propagateChange(-1, true);
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
	else if (rType == RefinementCrude)
	{
		Timer timer("positions to density", true);
		saveAtomPositions();
		double before = -scoreWithMap(ScoreTypeCorrel, target);
		int start = monomerBegin() + 0;
		int end = monomerEnd() - 0;

		for (int i = start; i < end; i += 3)
		{
			refineAroundMonomer(i, target);
		}
		double after = -scoreWithMap(ScoreTypeCorrel, target);

		std::cout << "Overall CC " << before * 100 << " to "
		<< after * 100 << "%." << std::endl;
		timer.report();

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
	_graphName = graphName;
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

		if (!c)
		{
			continue;
		}

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

void Polymer::refitBackbone(int start_, int end_)
{
	std::cout << "Refit backbone function" << std::endl;
	
	/* We can't straddle the anchor */
	if (start_ < _anchorNum && end_ > _anchorNum ||
	    end_ > _anchorNum && start_ < _anchorNum)
	{
		warn_user("We cannot straddle the anchor.");
		return;
	}

	downWeightResidues(start_, end_, 0);
	
	AtomPtr endAtom = getMonomer(end_)->findAtom("C");
	AtomPtr stAtom = getMonomer(start_)->findAtom("N");

	vec3 pos = stAtom->getAbsolutePosition();
	Options::getRuntimeOptions()->focusOnPosition(pos);
	
	/* Get new map for this */
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	Crystal::vsConcludeRefinement(&*ToParserPtr(crystal));
	
	/* Lower the number of samples */
	int samples = getAnchorModel()->getFinalPositions().size();
	Options::setNSamples(NULL, 0);
	refreshPositions();
	
	int start = start_;
	int end = end_;
	
	/* Swap if in wrong order */
	if ((start > _anchorNum && end < start)	||
	    (start < _anchorNum && end > start))
	{
		int tmp = start;
		start = end;
		end = tmp;
	}
	
	/* Split the bond */
	if (start < _anchorNum)
	{
		endAtom = getMonomer(end)->findAtom("N");
		stAtom = getMonomer(start)->findAtom("C");
	}
	
	BondPtr endBond = ToBondPtr(endAtom->getModel());
	BondPtr stBond = ToBondPtr(stAtom->getModel());
	
	endBond->setSplitBlock();
	BondPtr dupl = stBond->splitBond();

	Refitter fit(dupl, endBond, start > _anchorNum);
	fit.refit();

	downWeightResidues(start_, end_, 1);

	std::cout << "Return number of samples" << std::endl;
	/* Repair the number of samples */
	Options::setNSamples(NULL, samples);
	refreshPositions();
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

double Polymer::overfitTest(int round)
{
	int increment = (monomerCount() - beginMonomer()->first) / 10;
	int test = increment / 2 + beginMonomer()->first + round * increment;

	MonomerPtr monomer = getMonomer(test);
	std::cout << "Overfit test " << round << " on polymer "
	<< getChainID() << ", residue " << monomer->getResidueNum() << std::endl;
	
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	
	if (!monomer)
	{
		return std::nan(" ");
	}
	
	BackbonePtr bone = monomer->getBackbone();
	
	double st_obs = bone->scoreWithMap(ScoreTypeAddDensity, crystal,
	                                 false, MapScoreFlagNone);
	double st_diff = bone->scoreWithMap(ScoreTypeAddDensity, crystal,
	                                 false, MapScoreFlagDifference);
	
	double weight = 0;

	while (true)
	{
		for (int i = 0; i < bone->atomCount(); i++)
		{
			AtomPtr atom = bone->atom(i);
			atom->setWeightOnly(weight);
		}

		crystal->silentConcludeRefinement();

		double obs = bone->scoreWithMap(ScoreTypeAddDensity, crystal,
		                                false, MapScoreFlagNone);

		double diff = bone->scoreWithMap(ScoreTypeAddDensity, crystal,
		                                 false, MapScoreFlagDifference);

		double ratio = obs / st_obs;
		double diffratio = (diff - st_diff) / (st_obs - st_diff);
		std::cout << weight << ", " << ratio << ", " << diffratio
		<< ", " << ratio + diffratio << std::endl;
		
		if (ratio > 2)
		{
			break;
		}
		
		weight += 0.4;
	}

	for (int i = 0; i < bone->atomCount(); i++)
	{
		AtomPtr atom = bone->atom(i);
		atom->setWeightOnly(1);
	}

	return weight;
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
			
			if (bond->hasWhack())
			{
				WhackPtr whack = bond->getWhack();
				double k = Whack::getKick(&*whack);
				double w = Whack::getWhack(&*whack);
				k*= value;
				w *= value;
				Whack::setKick(&*whack, k);
				Whack::setWhack(&*whack, w);

			}
			else
			{
				double kick = Bond::getKick(&*bond);
				kick *= value;
				Bond::setKick(&*bond, kick);
			}
			
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

AnchorPtr Polymer::getAnchorModel()
{
	MonomerPtr anchoredRes = getMonomer(getAnchor());
	if (!anchoredRes)
	{
		return AnchorPtr();
	}

	ModelPtr model = anchoredRes->findAtom("N")->getModel();

	return ToAnchorPtr(model);
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
	std::cout << "\tB factor (Å^2): " << std::setprecision(4) <<
	bSum << std::endl;
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
	int flexCountBackbone = 0;
	
	int torsionCountSidechain = 0;
	int angleCountSidechain = 0;
	int flexCountSidechain = 0;

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
			BondPtr bond = ToBondPtr(a->getModel());
			
			if (bond->isTorsionRefinable() && bond->getRefineFlexibility())
			{
				int add = 1;
				
				if (bond->hasWhack())
				{
					add = 2;
				}

				if (back)
				{
					flexCountBackbone += add;
				}
				else
				{
					flexCountSidechain += add;
				}
			}

			if (bond->isTorsionRefinable())
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

			if (bond->getRefineBondAngle())
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
	
	int flexChain = 6 * 2;
	
	_positionalParams = torsionCountBackbone + torsionCountSidechain 
	+ angleCountBackbone + angleCountSidechain;
	
	_flexibilityParams = flexCountBackbone + flexCountSidechain + flexChain;
	
	std::cout << std::endl;
	std::cout << "|----------------" << std::endl;
	std::cout << "| Bond-based parameter groups (" << getChainID() << "):" << std::endl;
	std::cout << "|----------------" << std::endl;
	std::cout << "|          |   Backbone |  Sidechain |" << std::endl;
	std::cout << "| Torsion  |  " << std::setw(9) << torsionCountBackbone;
	std::cout << " |  " << std::setw(9) << torsionCountSidechain;
	std::cout << " | (positional)" << std::endl;
	std::cout << "| Angle    |  " << std::setw(9) << angleCountBackbone;
	std::cout << " |  " << std::setw(9) << angleCountSidechain;
	std::cout << " | (positional)" << std::endl;
	std::cout << "| Kicks    |  " << std::setw(9) << flexCountBackbone;
	std::cout << " |  " << std::setw(9) << flexCountSidechain;
	std::cout << " | (flexibility)" << std::endl;
	std::cout << "| Chain    |  " << std::setw(9) << flexChain;
	std::cout << " |  " << std::setw(9) << " ";
	std::cout << " | (whole molecule flex)" << std::endl;
	std::cout << "|----------------" << std::endl;
	std::cout << std::endl;
}

void Polymer::refineGlobalFlexibility()
{
	std::cout << "Optimising anchor shifts to match the electron density." << std::endl;

	Timer timer("anchor fit", true);
	
	AnchorPtr anchor = getAnchorModel();
	
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

