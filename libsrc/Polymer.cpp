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
#include "KeyPoints.h"
#include "Fibonacci.h"
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
	
	if (_startB < 0)
	{
		_startB = 20;
	}
	
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

bool Polymer::refineLocalFlexibility(bool magic)
{
	if (_anchorNum == -INT_MAX)
	{
		return false;
	}

	Timer timer("flexibility refinement", true);
	FlexLocal local;
	local.refineMagic(magic);
	whack();
	local.setPolymer(shared_from_this(), _kickShift);
	local.refine();
	_kickShift = local.getShift();
	timer.report();
	
	bool ch = local.didChange();
	
	if (!_keyPoints)
	{
		return ch;
	}
	
	if (ch)
	{
		return ch;
	}

	ch |= _keyPoints->refineKeyPoints();

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
	
	AtomList atoms = mon->findAtoms("CA");
	AnchorPtr anchor = getAnchorModel();

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
	
	if (bond->hasWhack())
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
}

void Polymer::whack()
{
	std::cout << "Whacking chain " << getChainID() << std::endl;
	
	for (int i = getAnchor(); i < monomerEnd(); i++)
	{
		whackMonomer(getMonomer(i));
	}

	for (int i = getAnchor() - 1; i >= monomerBegin(); i--)
	{
		whackMonomer(getMonomer(i));
	}

	getAnchorModel()->forceRefresh();
	refreshPositions();
}

void Polymer::tieAtomsUp()
{
	if (_anchorNum == -INT_MAX)
	{
		return;
	}
	
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

	bool newly_tied = false;
	if (nModel->isAbsolute())
	{
		AnchorPtr newAnchor = AnchorPtr(new Anchor(ToAbsolutePtr(nModel)));
		newAnchor->setBFactor(_startB);
		newAnchor->setNeighbouringAtoms(prev_ca, prev_c, ca, c);
		n->setModel(newAnchor);
		newly_tied = true;
	}

	for (int i = _anchorNum; i < monomerEnd(); i++)
	{
		if (!newly_tied && i == _anchorNum)
		{
			continue;
		}

		if (getMonomer(i))
		{
			getMonomer(i)->tieAtomsUp();
		}
	}

	for (int i = _anchorNum - 1; i >= monomerBegin(); i--)
	{
		if (newly_tied && i == _anchorNum)
		{
			continue;
		}

		if (getMonomer(i))
		{
			getMonomer(i)->tieAtomsUp();
		}
	}
	
	BondPtr n2ca = ToBondPtr(ca->getModel());
	BondPtr ca2c = ToBondPtr(c->getModel());
	BondPtr n2c = ToBondPtr(prev_c->getModel());
	BondPtr c2ca = ToBondPtr(prev_ca->getModel());

	if (newly_tied)
	{
		ca2c->setHeavyAlign(prev_c);
		n2ca->setHeavyAlign(prev_ca);
		n2c->setHeavyAlign(c);
		c2ca->setHeavyAlign(ca);
	}
	
	n2ca->checkForSplits(shared_from_this());
	n2c->checkForSplits(shared_from_this());
	
	/* test */
	if (Options::getBondAngles() >= 4)
	{
		for (int i = 0; i < atomCount(); i++)
		{
			if (atom(i)->getModel()->isBond())
			{
				ToBondPtr(atom(i)->getModel())->setRefineBondAngle(true);
			}
		}
	}
	
	whack();
	setupKeyPoints();
	setTied();
}

void Polymer::setupKeyPoints()
{
	KeyPointsPtr kp = KeyPointsPtr(new KeyPoints());
	kp->setPolymer(shared_from_this());
	_keyPoints = kp;
}

void Polymer::removeAtom(AtomPtr atom)
{
	MonomerPtr mon = atom->getMonomer();
	
	if (mon)
	{
		mon->removeAtom(atom);
	}
	
	if (atom->getModel()->isBond())
	{
		BondPtr bond = ToBondPtr(atom->getModel());

		if (bond->hasWhack() && getAnchorModel())
		{
			getAnchorModel()->removeWhack(bond->getWhack());
		}
	}
	
	AtomGroup::removeAtom(atom);
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
	
	addParamType(ParamOptionMaxTries, 1.0);
	addParamType(ParamOptionTorsion, 0.5);
	addParamType(ParamOptionNumBonds, 12);
	addParamType(ParamOptionExtraAtoms, 2);
	addParamType(ParamOptionSVD, 1);

	refineToEnd(position, crystal, rType);
}

void Polymer::refineBackbone()
{
	int anchor = getAnchor();
	clearParams();

	refineBackboneFrom(anchor - 3);
	refineBackboneFrom(anchor + 3);
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
	
	double step = 1.0;

	std::cout << "Refining region " << getChainID() << "   " << 
	std::right << std::setw(3) << coreStart << " - " << std::setw(3) <<
	coreEnd << " |   " << std::flush;
	
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
	vec3 c = coreRegion->centroid();
	eModel->setRotCentre(c);
//	refineShiftModel(target, eModel, coreRegion);

	coreRegion->saveScore();
	coreRegion->addParamType(ParamOptionNumBonds, (diff + 1) * 3);
	coreRegion->addParamType(ParamOptionTorsion, step);
	coreRegion->addParamType(ParamOptionTwist, -step);
//	coreRegion->addParamType(ParamOptionShift, 1);
	coreRegion->addParamType(ParamOptionMaxTries, 1);

	coreRegion->refine(target, RefinementCrude);
	
	/*
	coreRegion->clearParams();
	coreRegion->addParamType(ParamOptionNumBonds, (diff + 1) * 3);
	coreRegion->addParamType(ParamOptionMaxTries, 1);
	coreRegion->refine(target, RefinementCrude);
	*/

	coreRegion->saveAtomPositions();
	eModel->clearTwists();
	
	bool coversAnchor = (anchor >= coreStart && anchor <= coreEnd);
	
	int nTerm = -1;
	int nStart = anchor - 1;
	int cStart = anchor + 1;
	int cTerm = -1;
	
	int pre_step = 3;
	int post_step = 12;
	
	if (coreStart < anchor && !coversAnchor)
	{
		nStart = coreEnd + pre_step;
		if (nStart >= anchor)
		{
			nStart = anchor - 1;
		}
		
		nTerm = coreStart - post_step;
	}	
	
	if (coreEnd > anchor && !coversAnchor)
	{
		cStart = coreStart - pre_step;
		if (cStart > anchor)
		{
			cStart = coreEnd;
		}
		
		cTerm = coreEnd + post_step;
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

	double ccdiff = coreRegion->getImprovement();
	ccdiff *= 100;
	print_cc_diff(ccdiff, -1);
	std::cout << std::endl;

	return;
}

void Polymer::refineToEnd(int monNum, CrystalPtr target, RefinementType rType)
{
	int start = monNum;
	int end = (monNum < _anchorNum) ? monomerBegin() - 1: monomerEnd() + 1;
	
	refineRange(start, end, target, rType);
}

AtomGroupPtr Polymer::monomerRange(int start, int end, bool side)
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
		
		if (side)
		{
			all->addAtomsFrom(monomer);
			continue;
		}

		BackbonePtr bone = monomer->getBackbone();

		if (!bone)
		{
			continue;
		}
		
		all->addAtomsFrom(bone);
	}
	
	return all;
}

double Polymer::refineRange(int start, int end, CrystalPtr target, 
                            RefinementType rType)
{
	int skip = (start < _anchorNum) ? -1 : 1;
	if (((_anchorNum > start && _anchorNum < end) ||
	    (_anchorNum > end && _anchorNum < start)))
	{
		return 0;
		shout_at_helen("Trying to refine a range which straddles the anchor.");	
	}

	if (end > start && skip < 0)
	{
		int tmp = end;
		end = start;
		start = tmp - 1;
	}
	else
	{
		end++;
	}
	
	Timer timer("refine range", true);

	int count = 0;
	double startCCAve = 0;
	std::map<MonomerPtr, vec3> preScores;

	std::cout << "Refining chain " << getChainID();
	std::cout  << " from residue " << start << " towards ";
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
		
		if (rType == RefinementFine && getMonomer(i + skip))
		{
			BackbonePtr second = getMonomer(i + skip)->getBackbone();
			vec3 centre = bone->centroid();
			Options::getRuntimeOptions()->focusOnPosition(centre, 24);
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
		
		print_cc_diff(backdiff, 20);
		std::cout << "| " << std::flush;
		print_cc_diff(sidediff, -1);
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

void Polymer::refineShiftModel(CrystalPtr target, ExplicitModelPtr e,
                               AtomGroupPtr atoms)
{
	if (e->isAnchor())
	{
		return;
	}

	setupNelderMead();
	setCrystal(target);
	setCycles(100);
	setScoreType(ScoreTypeCorrel);

	setSilent(true);
	addParamType(ParamOptionShift, 1);
	addSampledAtoms(atoms);

	addTwistShift(e, atoms);

	sample();
}

void Polymer::refineAnchorPosition(CrystalPtr target)
{
	bool changed = true;
	int count = 0;
	saveAtomPositions();
	
	AtomGroupPtr surround = monomerRange(_anchorNum - 2, _anchorNum + 2,
	                                     false);
	
	AtomPtr n = getAnchorModel()->getNAtom();
	AtomPtr c = getAnchorModel()->getCAtom();
	
	BondPtr n_next = ToBondPtr(n->getModel());
	BondPtr c_next = ToBondPtr(c->getModel());
	
	while (changed && count < 1)
	{
		count++;
		setupNelderMead();
		setCrystal(target);
		setCycles(200);
		setScoreType(ScoreTypeCorrel);

		setSilent(true);
		addParamType(ParamOptionNumBonds, 8);
		
		setupThoroughSet(n_next, false);
		setupThoroughSet(c_next, false);
		addAnchorParams(getAnchorModel());

		changed = sample();
	}

	surround->saveAtomPositions();
}

void Polymer::refine(CrystalPtr target, RefinementType rType)
{
	if (_anchorNum == -INT_MAX)
	{
		return;
	}
	
	
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
		refineAnchorPosition(target);
		
		int step = 3;

		for (int i = getAnchor() - 1; i >= start; i -= step)
		{
			refineFromFarAroundMonomer(i, target);
		}
		
		for (int i = getAnchor(); i <= end; i += step)
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
	std::cout << "Refining chain " << getChainID() << 
	" from anchor to N-terminus... (" << monomerBegin() << ")";
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
	std::cout << "Refining chain " << getChainID() << 
	" from anchor to C-terminus... (" << monomerEnd() << ")";
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
	CSVPtr sidechainCsv = CSVPtr(new CSV(3, "resnum", "oldB", "newB"));
	CSVPtr ccCsv = CSVPtr(new CSV(3, "resnum", "cc", "rf"));
	CrystalPtr target = Options::getActiveCrystal();

	bool extra = Options::makeDiagnostics();

	for (int i = monomerBegin(); i < monomerEnd(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}

		BackbonePtr backbone = getMonomer(i)->getBackbone();
		SidechainPtr sidechain = getMonomer(i)->getSidechain();
		double cc = 0; double rf = 0;

		if (extra)
		{
			cc = -sidechain->scoreWithMap(ScoreTypeCorrel, target);
			rf = sidechain->scoreWithMap(ScoreTypeRFactor, target);
		}

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
		ccCsv->addEntry(3, value, cc, rf);

		if (caModel->getClassName() == "Bond")
		{
			BondPtr caBond = boost::static_pointer_cast<Bond>(caModel);
			double meanSq = caBond->getMeanSquareDeviation();

			double posDisp = ca->posDisplacement();
			double sideDisp = sidechain->getAverageDisplacement();

			csv->addEntry(5, value, meanSq, ca->getInitialBFactor(),
			              posDisp, sideDisp);
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
	plotMap["yMax0"] = "100";
	plotMap["yMax1"] = "100";

	plotMap["xTitle0"] = "Residue number";
	plotMap["yTitle0"] = "B equivalent";
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
		plotMap["yTitle0"] = "B equivalent";
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
		plotMap["filename"] = "density_" + graphName;
		plotMap["height"] = "700";
		plotMap["width"] = "1200";
		plotMap["xHeader0"] = "resnum";
		plotMap["yHeader0"] = "cc";
		plotMap["colour0"] = "black";
		plotMap["yMin0"] = "0";
		plotMap["yMax0"] = "1";

		plotMap["xTitle0"] = "Sidechain number";
		plotMap["yTitle0"] = "Value of CC or RF";
		plotMap["style0"] = "line";

		ccCsv->setSubDirectory("density_fit");
		ccCsv->setSubDirectory("density_fit");
		ccCsv->plotPNG(plotMap);
		ccCsv->writeToFile("density_fit" + graphName + ".csv");

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
	CSVPtr kicks = CSVPtr(new CSV(3, "res", "whack", "kick"));

	for (int i = monomerBegin(); i < monomerEnd(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}
		
		bool forwards = (i > getAnchor());

		AtomList phiAtoms = getMonomer(i)->findAtoms("N");
		AtomList psiAtoms = getMonomer(i)->findAtoms("C");
		AtomList caAtoms = getMonomer(i)->findAtoms("CA");
		
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
			
			csv->addEntry(3, (double)i, rad2deg(tPhi), rad2deg(tPsi));
		}
		
		for (int j = 0; j < caAtoms.size() && j < 1; j++)
		{
			ModelPtr ca = caAtoms[j]->getModel();
			if (!ca->isBond() || !ToBondPtr(ca)->hasWhack())
			{
				continue;
			}
			
			WhackPtr whack = ToBondPtr(ca)->getWhack();
			double k = Whack::getKick(&*whack);
			double w = Whack::getWhack(&*whack);
			
			kicks->addEntry(3, (double)i, w, k); 
		}
	}

	csv->writeToFile("ramachandran.csv");
	kicks->writeToFile("chainflex_" + getChainID() + ".csv");
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
	Options::getRuntimeOptions()->focusOnPosition(pos, 24);
	
	/* Get new map for this */
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	crystal->wrapUpRefinement();
	
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

	int anchorRes = -INT_MAX;
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
		
		if (getMonomer(i - 1) && getMonomer(i - 1)->conformerCount() > 1)
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

	if (anchorRes == -INT_MAX)
	{
		std::cout << "Not anchoring " << getChainID() << std::endl;
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
	std::cout << "\tDerived B equivalent (Å^2): " << std::setprecision(4) <<
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
	
	int flexChain = 9;
	
	if (getAnchorModel())
	{
		flexChain += 15;
	}

	
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
		
		AtomList cbs = getMonomer(i)->findAtoms("CB");
		
		for (int j = 0; j < cbs.size(); j++)
		{
			_allBackbones->addAtom(cbs[j]);
		}
	}
	
	_allBackbones->setName(_name + "_backbone");
	return _allBackbones;
}

void Polymer::attachTargetToRefinement(RefinementStrategyPtr strategy,
                                       FlexGlobal &target, bool isotropy)
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	AtomGroupPtr allBackbone = getAllBackbone();
	target.setAtomGroup(allBackbone);
	target.setCrystal(crystal);
	
	strategy->setEvaluationFunction(FlexGlobal::score, &target);
	FlexGlobal::score(&target);
}

void Polymer::addProperties()
{
	Molecule::addProperties();

	addIntProperty("anchor_res", &_anchorNum);
	addChild("keypoints", _keyPoints);

	for (int i = monomerBegin(); i < monomerEnd(); i++)
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

	if (category == "keypoints")
	{
		KeyPointsPtr kp = ToKeyGroupPtr(object);
		_keyPoints = kp;
	}

	Molecule::addObject(object, category);
}

void Polymer::postParseTidy()
{
	Molecule::postParseTidy();
	//applyTranslationTensor();
	
	_name = "polymer_" + getChainID();
}

bool Polymer::hasResidue(int resNum)
{
	return (resNum >= monomerBegin() && resNum <= monomerEnd());
}

void Polymer::refineMotions()
{
	AnchorPtr anch = getAnchorModel();
	
	for (int i = 0; i < anch->motionCount(); i++)
	{
		MotionPtr mot = anch->getMotion(i);
		mot->refine();
	}
}

void Polymer::removeIntramolecularMotion()
{
	std::cout << "Wiping intramolecular motions for Polymer " << getChainID() 
	<< std::endl;
	for (int i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getModel()->isBond())
		{
			BondPtr b = ToBondPtr(atom(i)->getModel());
			Bond::setKick(&*b, 0.);
			
			if (b->hasWhack())
			{
				WhackPtr w = b->getWhack();
				Whack::setKick(&*w, 0.);
				Whack::setWhack(&*w, 0.);
			}
		}
	}

	getAnchorModel()->forceRefresh();
	refreshPositions();
}

void Polymer::redefineMotion()
{
	_allBackbones = AtomGroupPtr();
	AnchorPtr anch = getAnchorModel();
	
	for (int i = 0; i < anch->motionCount(); i++)
	{
		MotionPtr mot = anch->getMotion(i);
		mot->updateAtoms();
	}
	
	anch->forceRefresh();
	refreshPositions();
}

void Polymer::resetMotion()
{
	redefineMotion();

	AnchorPtr anch = getAnchorModel();
	
	for (int i = 0; i < anch->motionCount(); i++)
	{
		MotionPtr mot = anch->getMotion(i);
		mot->reset();
	}

	refreshPositions();
}

void Polymer::resetSidechains()
{
	for (int i = monomerBegin(); i < monomerEnd(); i++)
	{
		MonomerPtr mon = getMonomer(i);
		if (!mon)
		{
			continue;
		}
		SidechainPtr side = mon->getSidechain();
		side->reset();
	}

	refreshPositions();
}
