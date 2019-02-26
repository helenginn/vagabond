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
#include "Gradiator.h"
#include "Anchor.h"
#include "Whack.h"
#include "RefinementGridSearch.h"
#include "RefinementNelderMead.h"
#include "ParamBand.h"
#include "FlexGlobal.h"
#include "Reflex.h"
#include "Timer.h"
#include <map>
#include <iomanip>

FlexLocal::FlexLocal()
{
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

	_svd = new SVDBond(_bondEffects, _bonds, _atoms);
	_svd->setPolymer(_polymer);
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
	nelder->setCycles(200);
	nelder->setVerbose(true);	
	nelder->setSilent(true);

	nelder->setEvaluationFunction(getScore, this);

	_svd->addToStrategy(nelder, _negMult);

	nelder->refine();

	double val = (1 + getScore(this)) * 100.;
	val = nelder->improvement();

	if (nelder->didChange())
	{
		std::cout << std::setw(3) << val << 
		"% improved. ... done. ";
	}
	else
	{
		std::cout << " not improved.   ... done. ";
	}

	timer.quickReport();
	std::cout << std::endl;
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
		bondTest();
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
		local->_workspace.selectAtoms = bb->getAtoms();
		local->_prepared = true;
		AtomGroup::scoreWithMapGeneral(&local->_workspace);
	}
	
	if (local->_svd)
	{
		local->_svd->applyParameters();
	}

	local->_polymer->propagateChange();
	
	double score = AtomGroup::scoreWithMapGeneral(&local->_workspace);
	return score;
	
}

void FlexLocal::bondTest()
{
	Gradiator *g = _svd->getGradiator();
	SVDBond *svd = _svd;
	_svd = NULL;
	
	if (!g) return;
	
	double baseCC = getScore(this);
	CSVPtr csv = CSVPtr(new CSV(3, "bond", "whack", "kick"));
	
	for (size_t i = 0; i < g->whackCount(); i++)
	{
		WhackPtr w = g->getWhack(i);
		double val = Whack::getWhack(&*w);
		Whack::setWhack(&*w, val + 0.0001);
		double whackCC = getScore(this);
		Whack::setWhack(&*w, val);
		val = Whack::getKick(&*w);
		Whack::setKick(&*w, val + 0.0001);
		double kickCC = getScore(this);
		Whack::setKick(&*w, val);
		
		double wDiff = whackCC - baseCC;
		double kDiff = kickCC - baseCC;
		csv->addEntry(3, (double)i, wDiff, kDiff);
	}
	
	csv->setSubDirectory("local_flex");
	csv->writeToFile("bond_test.csv");
	
	_svd = svd;
}
