// Vagabond
// Copyright (C) 2017-2018 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

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
#include "Anisotropicator.h"
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
#include "RefinementLBFGS.h"
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

	for (int i = monomerBegin(); i < monomerEnd(); i++)
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

bool Polymer::refineLocalFlexibility()
{
	Timer timer("flexibility refinement", true);
	FlexLocal local;
	whack();
	local.setPolymer(shared_from_this(), _kickShift);
	local.refine();
	_kickShift = local.getShift();
	timer.report();
	
	bool ch = local.didChange();
	return ch;
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

		if (next->getAtom()->getResidueNum() != num)
		{
			break;
		}

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
	
	if (monomerCount() == 1)
	{
		warn_user("Not tying up " + getChainID() + " - only one monomer");
		return;
	}

	checkChainContinuity();

	AtomPtr n = getMonomer(_anchorNum)->findAtom("N");
	
	if (!getMonomer(_anchorNum - 1))
	{
		warn_user("Not tying up " + getChainID() + " - cannot find atoms "\
		"surrounding anchor point.");
		return;
	}
	
	/* Specify heavy alignment atoms around the anchor point */
	AtomPtr prev_c = getMonomer(_anchorNum - 1)->findAtom("C");
	AtomPtr prev_ca = getMonomer(_anchorNum - 1)->findAtom("CA");
	AtomPtr ca = getMonomer(_anchorNum)->findAtom("CA");
	AtomPtr c = getMonomer(_anchorNum)->findAtom("C");

	ModelPtr nModel = n->getModel();

	if (nModel->isAbsolute())
	{
		AnchorPtr newAnchor = AnchorPtr(new Anchor(ToAbsolutePtr(nModel)));
		newAnchor->setBFactor(_startB);
		newAnchor->setNeighbouringAtoms(prev_ca, prev_c, ca, c);
		n->setModel(newAnchor);
	}

	for (int i = _anchorNum; i < monomerEnd(); i++)
	{
		if (getMonomer(i))
		{
			getMonomer(i)->tieAtomsUp();
		}
	}

	for (int i = _anchorNum - 1; i >= monomerBegin(); i--)
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
	for (int i = monomerBegin(); i < monomerEnd(); i++)
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

void refineRegion(AtomGroupPtr region, CrystalPtr target, double light)
{
	region->addParamType(ParamOptionTorsion, light);
	region->addParamType(ParamOptionBondAngle, light);
	region->addParamType(ParamOptionNumBonds, 4);
	region->refine(target, RefinementSavedPos);
}

void Polymer::refineFromFarAroundMonomer(int central, CrystalPtr target)
{
	MonomerPtr monomer = getMonomer(central);
	if (!monomer)
	{
		return;
	}
	
	int pad = 2;
	int coreStart = central - pad;
	int coreEnd = central + pad;
	
	refineFromFarRegion(coreStart, coreEnd, target);
}
	
	
void Polymer::refineFromFarRegion(int coreStart, int coreEnd,
                                  CrystalPtr target)
{
	size_t diff = coreEnd - coreStart;
	
	AtomGroupPtr coreRegion = monomerRange(coreStart, coreEnd);
	
	if (!coreRegion || coreRegion->atomCount() == 0)
	{
		return;
	}
	
	int anchor = getAnchor();
	
	double step = 0.5;

	std::cout << "Refining core region, " << getChainID() << " " << coreStart
	<< " to " << coreEnd << std::flush;
	
	AtomList top = coreRegion->topLevelAtoms();
	
	if (!top.size())
	{
		return;
	}
	
	ModelPtr model = top[0]->getModel();
	
	ExplicitModelPtr eModel;
	if (model->isBond())
	{
		eModel = ToBondPtr(model)->getParentModel();
	}
	else if (model->isAnchor())
	{
		eModel = ToAnchorPtr(model);
	}
	else
	{
		return;
	}

	coreRegion->makeBackboneTwists(eModel);

	coreRegion->addParamType(ParamOptionNumBonds, (diff + 1) * 3);
	coreRegion->addParamType(ParamOptionTorsion, step);
	coreRegion->addParamType(ParamOptionTwist, step);
	coreRegion->addParamType(ParamOptionBondAngle, step);
	coreRegion->addParamType(ParamOptionMaxTries, 1);
	coreRegion->addIncludeForRefinement(coreRegion);
	coreRegion->refine(target, RefinementCrude);

	coreRegion->saveAtomPositions();
	eModel->clearTwists();
	
	bool coversAnchor = (anchor >= coreStart && anchor <= coreEnd);
	
	if (coversAnchor)
	{
		refineAnchorPosition(target);
	}
	
	std::cout << "." << std::flush;
	
	int nTerm = -1;
	int nStart = anchor - 1;
	int cStart = anchor + 1;
	int cTerm = -1;
	
	if (coreStart < anchor && !coversAnchor)
	{
		nStart = coreEnd + 3;
		if (nStart >= anchor)
		{
			nStart = anchor - 1;
		}
		
		nTerm = coreStart - 10;
	}	
	
	if (coreEnd > anchor && !coversAnchor)
	{
		cStart = coreStart - 3;
		if (cStart > anchor)
		{
			cStart = coreEnd;
		}
		
		cTerm = coreEnd + 10;
	}
	
	AtomGroupPtr leftRegion = monomerRange(nTerm, nStart);
	AtomGroupPtr rightRegion = monomerRange(cStart, cTerm);

	if (coreEnd < anchor || coversAnchor)
	{
		refineRegion(leftRegion, target, step);
	}

	if (coreStart >= anchor || coversAnchor)
	{
		refineRegion(rightRegion, target, step);
	}

	std::cout << " Displacement by " << coreRegion->getAverageDisplacement()
	<< " Å." << std::endl;

}

void Polymer::refineToEnd(int monNum, CrystalPtr target, RefinementType rType)
{
	int start = monNum;
	int end = (monNum < _anchorNum) ? monomerBegin() : monomerEnd();
	
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

	getAnchorModel()->propagateChange(20, true);
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

		for (int i = getAnchor() - 1; i >= start; i -= 3)
		{
			refineFromFarAroundMonomer(i, target);
		}
		
		for (int i = getAnchor(); i <= end; i += 3)
		{
			refineFromFarAroundMonomer(i, target);
		}

		getAnchorModel()->propagateChange(-1, true);

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

	for (int i = _anchorNum - 1; i >= monomerBegin(); i--)
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

	for (int i = _anchorNum; i < monomerEnd(); i++)
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

	for (int i = monomerBegin(); i < monomerEnd(); i++)
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

	for (int i = monomerBegin(); i < monomerEnd(); i++)
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
	for (int i = monomerBegin(); i < monomerEnd(); i++)
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
	double count = 0;
	vec3 sum = centroid();

	int anchorRes = -1;
	double lowestLength = FLT_MAX;

	for (int i = monomerBegin() + 1; i < monomerEnd() - 1; i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}
		
		if (getMonomer(i)->conformerCount() > 1)
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

	for (int i = monomerBegin(); i < monomerEnd(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}	
		
		hydrogenator.setMonomer(getMonomer(i));
		hydrogenator.hydrogenate();
	}
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
				else if (isWhacking() && 
				         (a->isBackbone() || a->isBackboneAndSidechain()))
				{
					add = 0;
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
	std::cout << "Optimising whole molecule movements, chain " << 
	getChainID() << "." << std::endl;

	Timer timer("anchor fit", true);
	
	AnchorPtr anchor = getAnchorModel();
	
	{
		FlexGlobal target;
//		RefinementLBFGSPtr lbfgs = RefinementLBFGSPtr(new RefinementLBFGS());
		NelderMeadPtr lbfgs = NelderMeadPtr(new RefinementNelderMead());
		attachTargetToRefinement(lbfgs, target);
		lbfgs->setJobName("translation");
		lbfgs->setCycles(24);

		anchor->addTranslationParameters(lbfgs);
		lbfgs->refine();
	}

	{
		FlexGlobal target;
		NelderMeadPtr lbfgs = NelderMeadPtr(new RefinementNelderMead());
		lbfgs->setCycles(24);
		lbfgs->setJobName("libscrew");
		attachTargetToRefinement(lbfgs, target);

		anchor->addLibrationParameters(lbfgs);
//		anchor->addScrewParameters(lbfgs);
		lbfgs->refine();
	}

	timer.report();
}

AtomGroupPtr Polymer::getAllBackbone()
{
	if (_allBackbones)
	{
		return _allBackbones;
	}
	
	_allBackbones = AtomGroupPtr(new AtomGroup());

	for (int i = monomerBegin(); i < monomerEnd(); i++)
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

void Polymer::addProperties()
{
	Molecule::addProperties();

	addIntProperty("anchor_res", &_anchorNum);

	for (int i = monomerBegin(); i < monomerEnd(); i++)
	{
		std::cout << i << std::endl;
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
	//applyTranslationTensor();
	
}

mat3x3 Polymer::fitEllipsoid()
{
	std::vector<vec3> points;

	for (int i = monomerBegin(); i < monomerEnd(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}
		
		AtomPtr ca = getMonomer(i)->findAtom("CA");
		points.push_back(ca->getAbsolutePosition());
	}

	Anisotropicator trop;
	trop.setPoints(points);
	
	mat3x3 basis = trop.getTensor();
	getAnchorModel()->setPolymerBasis(basis);

	return basis;
}

