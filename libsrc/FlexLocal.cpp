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
	_svd->performSVD();
	
	std::cout << _svd->numClusters() << " clusters.";
	std::cout << "              ... done. ";
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

		/* Testing the SVD stuff */
		svd();
		
		for (int j = 0; j < 1; j++)
		{
			std::random_shuffle(_paramBands.begin(), _paramBands.end());

			std::cout << "| " << j + 1 << ". Refining bond clusters... " 
			<< std::flush;
			Timer timer;
			int limit = 5;

			NelderMeadPtr nelder = NelderMeadPtr(new RefinementNelderMead());
			nelder->setCycles(200);
			nelder->setVerbose(true);	

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
		
		std::cout << "---------------------------------------------------------"
		<< std::endl;

		clear();

		/* Flip the sign of the kicks next time */
		_negMult *= -1;		
		_run++;
	}
}

void FlexLocal::clear()
{
	_atomTargets.clear();
	_atomOriginal.clear();
	_bondEffects.clear();
	_atoms.clear();
	_bonds.clear();
	_reorderedBonds.clear();
	_b2bDiffs.clear();
	_degrees.clear();
	_bondClusterIds.clear();
	_bbCCs.clear();
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

double FlexLocal::bondRelationship(BondPtr bi, BondPtr bj)
{
	std::vector<double> xs, ys;
	
	for (int i = 0; i < _atoms.size(); i++)
	{
		AtomPtr a = _atoms[i];
		double flexi = _bondEffects[bi][a];
		double flexj = _bondEffects[bj][a];
		
		if (fabs(flexi) < 1e-6 || fabs(flexj) < 1e-6)
		{
			continue;
		}
		
		xs.push_back(flexi);
		ys.push_back(flexj);
	}
	
	if (xs.size() == 1)
	{
		return nan(" ");
	}

	double correl = correlation(xs, ys);
	
	return correl;
}

bool less_than(const BondDegree b1, const BondDegree b2)
{
	return (b2.degree > b1.degree);
}

void FlexLocal::createClustering()
{
	CSVPtr csv = CSVPtr(new CSV(3, "bond_i", "bond_j", "bondcc"));

	double sumCC = 0;
	double count = 0;

	Timer timer;
	std::cout << "| 2. Calculating bond-bond correlation pairs... " 
	<< std::flush;

	/* Once all the bond effects are calculated, 
	 * 	give them pair-wise correlations */

	for (int i = 0; i < _bonds.size(); i++)
	{
		BondPtr bi = _bonds[i];
		
		int degree = 0;

		for (int j = 0; j < _bonds.size(); j++)
		{
			BondPtr bj = _bonds[j]; // giggle
			double cc = 0;

			cc = bondRelationship(bi, bj);
			
			if (cc != cc)
			{
				cc = 0;
			}
			
			_bbCCs[bi][bj] = cc;
			_bbCCs[bj][bi] = cc;

			sumCC += cc;
			count++;
			degree += (cc > _threshold);
			double reach = (cc > _threshold);
			
			csv->addEntry(3, (double)i, (double)j, cc);
		}
		
		BondDegree bd;
		bd.index = i;
		bd.degree = degree;
		_degrees.push_back(bd);
	}
	
	std::cout << "... done. ";
	timer.quickReport();
	std::cout << std::endl;

	csv->setSubDirectory("local_flex");
	csv->writeToFile(_polymer->getChainID() + "bond_matrix.csv");
	
	std::sort(_degrees.begin(), _degrees.end(), less_than);

}

double FlexLocal::targetForAtom(AtomPtr a)
{
	double target = _atomTargets[a];

	double transformed = target / _increment;
	
	return transformed;
}

double FlexLocal::actualAtomChange(AtomPtr a)
{
	double original = _atomOriginal[a];
	double current = a->getBFactor();

	return current - original;
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

void FlexLocal::reflex()
{
	if (!_useTarget)
	{
		return;
	}
	
	// this updates the atom B targets.
	Timer timer;
	std::cout << "| 0. Determining atom-flexibility targets... " << std::flush;
	Reflex reflex;
	reflex.setPolymer(_polymer);
	reflex.setPieceCount(1);
	reflex.calculate();
	std::cout << "   ... done. ";
	timer.quickReport();
	std::cout << std::endl;
}

void FlexLocal::scanBondParams()
{
	ExplicitModelPtr model = _polymer->getAnchorModel();
	int samples = model->getFinalPositions().size();

	if (samples > 40)
	{
		Options::setNSamples(NULL, 40);
	}

	_polymer->refreshPositions();

	Timer timer;
	std::cout << "| 1. Determining atom-flexibility effects... " << std::flush;
	
	CSVPtr atomtarg = CSVPtr(new CSV(3, "atom", "target", "increment"));
	
	for (int i = 0; i < _atoms.size(); i++)
	{
		double num = _atoms[i]->getResidueNum();
		double increment = targetForAtom(_atoms[i]);
		
		double target = _atoms[i]->getInitialBFactor();
		
		if (_useTarget)
		{
			target = _atoms[i]->getTargetB();
		}

		atomtarg->addEntry(3, num, target, increment);
	}

	atomtarg->setSubDirectory("local_flex");
	atomtarg->writeToFile("atom_targets_" + _polymer->getGraphName() + ".csv");

	{
		std::map<std::string, std::string> plotMap;
		plotMap["height"] = "700";
		plotMap["width"] = "1200";
		plotMap["xHeader0"] = "atom";
		plotMap["colour0"] = "black";
		plotMap["xTitle0"] = "Residue number";
		plotMap["yTitle0"] = "B factor target";
		plotMap["style0"] = "line";
		plotMap["yMin0"] = "-1";
		plotMap["yMax0"] = "2";
		plotMap["filename"] = "atom_increments_" + _polymer->getGraphName();
		plotMap["yHeader0"] = "increment";
		atomtarg->plotPNG(plotMap);
	}
	
	CSVPtr csv = CSVPtr(new CSV(3, "atom", "kick", "response"));
	int count = 0;
	
	std::map<int, double> sumsPerAtom;

	for (int i = 0; i < _bonds.size(); i++)
	{
		BondPtr b = _bonds[i];
		WhackPtr w = _bonds[i]->getWhack();
		
		if (!w)
		{
			continue;
		}

		double ki = Whack::getKick(&*w);
		double wh = Whack::getWhack(&*w);

		double add = _shift * _negMult;

		setBondParam(b, wh + add, ki + add);

		for (int j = 0; j < _atoms.size(); j++)
		{
			double nAtom = _atoms[j]->getResidueNum();
			double bnow = _atoms[j]->getBFactor();
			double bthen = _atomOriginal[_atoms[j]];

			double diff = bnow - bthen;

			if (_bondEffects.count(b) == 0)
			{
				_bondEffects[b] = AtomTarget();
			}

			_bondEffects[b][_atoms[j]] = diff;

			int intAtom = _atoms[j]->getResidueNum();

			if (sumsPerAtom.count(intAtom) == 0)
			{
				sumsPerAtom[intAtom] = 0;
			}

			sumsPerAtom[intAtom] += diff;

			csv->addEntry(3, (double)nAtom, (double)count, diff);

		}

		setBondParam(b, wh, ki);

		count++;
	}
	
	/* Output all sums per atom */

	CSVPtr sums = CSVPtr(new CSV(2, "resnum", "sum"));

	for (std::map<int, double>::iterator it = sumsPerAtom.begin();
	     it != sumsPerAtom.end(); it++)
	{
		sums->addEntry(2, (double)it->first, (double)it->second);
	}

	sums->setSubDirectory("local_flex");
	sums->writeToFile("sum_scan_" + _polymer->getGraphName() + ".csv");
	
	/* Output bond scan matrix */

	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = "bond_scan_" + _polymer->getGraphName()
	+ (_negMult > 0 ? "a" : "b");
	plotMap["height"] = "1000";
	plotMap["width"] = "1000";
	plotMap["xHeader0"] = "kick";
	plotMap["yHeader0"] = "atom";
	plotMap["zHeader0"] = "response";

	plotMap["xTitle0"] = "bond number";
	plotMap["yTitle0"] = "atom";

	plotMap["style0"] = "heatmap";
	plotMap["stride"] = i_to_str(_bonds.size());
	
	csv->setSubDirectory("local_flex");
	csv->writeToFile(plotMap["filename"] + ".csv");
	csv->plotPNG(plotMap);
	
	std::cout << "   ... done. ";
	timer.quickReport();
	std::cout << std::endl;

	Options::setNSamples(NULL, samples);
}

void FlexLocal::propagateWhack()
{
	ExplicitModelPtr anchor = _polymer->getAnchorModel();
	anchor->propagateChange(-1, true);
}

void FlexLocal::setBondParam(BondPtr b, double wh, double k)
{
	WhackPtr w = b->getWhack();
	if (!w) return;

	Whack::setWhack(&*w, wh);
	Whack::setKick(&*w, k);
	propagateWhack();
}

