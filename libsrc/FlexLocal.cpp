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
#include "Anchor.h"
#include "Whack.h"
#include <hcsrc/RefinementGridSearch.h>
#include <hcsrc/RefinementNelderMead.h>
#include "ParamBand.h"
#include "FlexGlobal.h"
#include <hcsrc/Timer.h>
#include <map>
#include <iomanip>

FlexLocal::FlexLocal()
{
	_magic = false;
	_prepared = false;
	_startB = 0;
	_shift = 0.05;
	_run = 0;
	_negMult = 1;
	_anchorB = 0;
	_threshold = 0.80;
	_increment = 5;
	_useTarget = true;
	_getter = Bond::getKick;
	_flexGlobal = NULL;
	_setter = Bond::setKick;
	_svd = NULL;
	_changed = false;
}

FlexLocal::~FlexLocal()
{
	if (_flexGlobal)
	{
		delete _flexGlobal;
		_flexGlobal = NULL;
	}
	
	if (_svd)
	{
		delete _svd;
		_svd = NULL;
	}
}

void FlexLocal::setPolymer(PolymerPtr pol, double shift)
{
	_polymer = pol;
	_shift = shift;
	_bb = _polymer->getAllBackbone();
}

void FlexLocal::svd()
{
	std::cout << "| 0. Performing SVD... " << std::flush;
	Timer timer;
	
	if (_svd)
	{
		delete _svd;
		_svd = NULL;
	}

	_svd = new SVDBond(_bonds, _atoms);
	_svd->performSVD();
	
	std::cout << _svd->numClusters() << " clusters.";
	std::cout << "              ... done. ";
	timer.quickReport();
	std::cout << std::endl;
}

void FlexLocal::refineClusters()
{
	std::cout << "| 1. Refining bond clusters... " 
	<< std::flush;
	Timer timer;

	NelderMeadPtr nelder = NelderMeadPtr(new RefinementNelderMead());
	nelder->setCycles(120);
	nelder->setVerbose(true);	
	nelder->setSilent(true);

	nelder->setEvaluationFunction(getScore, this);

	CrystalPtr crystal = Options::getActiveCrystal();

	_svd->addToStrategy(nelder, _negMult, _magic);
	nelder->refine();
	nelder->reportResult();

	timer.quickReport();
	std::cout << std::endl;
	_changed = nelder->changedSignificantly();
}

void FlexLocal::refine()
{
	for (int i = 0; i < 1; i++)
	{
		std::cout << "---------------------------------------------------------"
		<< std::endl;
		std::cout << "|  Refining flexibility for chain " 
		<< _polymer->getChainID() << " (cycle " << _run << ")";
		std::cout << std::endl;
		std::cout << "---------------------------------------------------------"
		<< std::endl;

		findAtomsAndBonds();
		svd();
		refineClusters();
		
		std::cout << "---------------------------------------------------------"
		<< std::endl;

		clear();
		_run++;
	}
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
	for (int i = 0; i < _polymer->atomCount(); i++)
	{
		AtomPtr a = _polymer->atom(i);

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
		AtomGroupPtr bb = local->_polymer->getAllBackbone();
		setup_space(&local->_workspace);
		local->_workspace.crystal = Options::getActiveCrystal();
		local->_workspace.selectAtoms = bb;
		local->_workspace.selectAtoms->addParamType(ParamOptionStep, 2);
		local->_prepared = true;
		AtomGroup::scoreWithMapGeneral(&local->_workspace);
	}
	
	if (local->_svd)
	{
		local->_svd->applyParameters();
	}

	local->_bb->refreshPositions();
	
	double score = AtomGroup::scoreWithMapGeneral(&local->_workspace);
	return score;
	
}

