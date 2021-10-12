//
//  FlexLocal.cpp
//  vagabond
//
//  Created by Helen Ginn, 2018
//  Copyright © 2018 Helen Ginn. All rights reserved.
//

#include "FlexLocal.h"
#include "CSV.h"
#include "Bond.h"
#include "Polymer.h"
#include <svdcmp.h>
#include "Monomer.h"
#include "SVDBond.h"
#include "Anchor.h"
#include <hcsrc/RefinementNelderMead.h>
#include "ParamBand.h"
#include <hcsrc/Timer.h>
#include <map>
#include <iomanip>

FlexLocal::FlexLocal()
{
	_prepared = false;
	_shift = 0.05;
	_run = 0;
	_svd = NULL;
	_changed = false;
	_torsionMode = false;
	_stream = &std::cout;
}

FlexLocal::~FlexLocal()
{
	if (_svd)
	{
		delete _svd;
		_svd = NULL;
	}
}

void FlexLocal::setPolymer(PolymerPtr pol)
{
	_polymer = pol;
	_stream = pol->getStream();
	_bb = _polymer->getAllBackbone();
}

void FlexLocal::svd()
{
	*_stream << "| 0. Performing SVD... " << std::flush;
	Timer timer;
	
	if (_svd)
	{
		delete _svd;
		_svd = NULL;
	}

	_svd = new SVDBond(_bonds, _atoms);
	_svd->setWriteMatrix(true);
	
	if (_torsionMode)
	{
		_svd->setDoTorsion(true);
	}

	_svd->performSVD();
	
	*_stream << _svd->numClusters() << " clusters.";
	*_stream << "              ... done. ";
	timer.quickReport();
	*_stream << std::endl;
}

void FlexLocal::refineClusters()
{
	*_stream << "| 1. Refining bond clusters... " << std::endl;
	Timer timer;

	NelderMeadPtr nelder = NelderMeadPtr(new RefinementNelderMead());
	nelder->setCycles(100);
	nelder->setVerbose(true);	
	nelder->setSilent(true);

	nelder->setEvaluationFunction(getScore, this);
	if (_polymer)
	{
		nelder->setStream(_polymer->getStream());
	}

	CrystalPtr crystal = Options::getActiveCrystal();

	_svd->addToStrategy(nelder, _shift, false);
	AnchorPtr anch = _polymer->getAnchorModel();
//	anch->addChainMultsToStrategy(nelder);
	
	nelder->refine();
	nelder->reportResult();

	timer.quickReport();
	*_stream << std::endl;
	_changed = nelder->changedSignificantly();
}

void FlexLocal::refine()
{
	for (int i = 0; i < 1; i++)
	{
		*_stream << "---------------------------------------------------------"
		<< std::endl;
		*_stream << "|  Refining flexibility for chain " 
		<< _polymer->getChainID() << " (cycle " << _run << ")";
		*_stream << std::endl;
		*_stream << "---------------------------------------------------------"
		<< std::endl;

		findAtomsAndBonds();
		svd();
		refineClusters();
		
		*_stream << "---------------------------------------------------------"
		<< std::endl;

		clear();
		_run++;
	}
}

void FlexLocal::attachToStrategy(RefinementStrategyPtr strategy,
                                 AtomGroupPtr group)
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	setCrystal(crystal);
	strategy->setStream(_stream);
	
	setAtomGroup(group);

	strategy->setVerbose(true);
	strategy->setCycles(60);

	strategy->setEvaluationFunction(FlexLocal::getScore, this);
	FlexLocal::getScore(this);

}

void FlexLocal::refineChainMults(AnchorPtr anch)
{
	Timer timer;
	*_stream << "| 1a. Refining kick spread... " 
	<< std::flush;
	NelderMeadPtr ref = NelderMeadPtr(new RefinementNelderMead());
	ref->setCycles(200);
	ref->setEvaluationFunction(getScore, this);
	ref->setVerbose(true);	
	ref->setSilent(true);
	anch->addChainMultsToStrategy(ref);

	ref->refine();
	ref->reportResult();


	timer.quickReport();
	*_stream << std::endl;
}

void FlexLocal::clear()
{
	_bondEffects.clear();
	_atoms.clear();
	_bonds.clear();
	_prepared = false;
}

void FlexLocal::findAtomsAndBonds()
{
	AtomGroupPtr chosen = _polymer;
	
	if (_bb)
	{
		chosen = _bb;
	}

	for (int i = 0; i < chosen->atomCount(); i++)
	{
		AtomPtr a = chosen->atom(i);

		ModelPtr m = a->getModel();
		
		if (!m->isBond())
		{
			continue;
		}
		
		BondPtr b = ToBondPtr(a->getModel());

		if (!b->getWhack())
		{
			continue;
		}
		
		if (a->getAtomName() != "CA")
		{
			continue;
		}

		_bonds.push_back(b);
		_atoms.push_back(a);
	}
}	

double FlexLocal::getScore(void *object)
{
	FlexLocal *local = static_cast<FlexLocal *>(object);

	if (!local->_prepared)
	{
		AtomGroupPtr bb = local->_bb;
		if (!bb || bb->atomCount() == 0)
		{
			bb = local->_polymer->getAllBackbone();
		}

		setup_space(&local->_workspace);
		local->_workspace.crystal = Options::getActiveCrystal();
		local->_workspace.selectAtoms = bb;
		local->_workspace.selectAtoms->addParamType(ParamOptionStep, 2);
		local->_prepared = true;
		AtomGroup::scoreWithMapGeneral(&local->_workspace);
	}

	local->_workspace.tBonds->setFine(true);
	local->_workspace.tMap->setFine(true);
	local->_workspace.tScore->setFine(true);
	
	if (local->_svd)
	{
		local->_svd->applyParameters();
	}

	local->_workspace.tBonds->start();
	local->_bb->refreshPositions();
	local->_workspace.tBonds->stop();
	
	double score = AtomGroup::scoreWithMapGeneral(&local->_workspace);
	return score;
}

void FlexLocal::recalculateConstant()
{
//	_workspace.recalc = true;
}

void FlexLocal::reportTimings()
{
	*_stream << std::endl;
	_workspace.tBonds->report(_stream);
	_workspace.tMap->report(_stream);
	_workspace.tScore->report(_stream);
}

