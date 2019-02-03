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
	_shift = 0.02;
	_run = 0;
	_window = 10;
	_direct = 0;
	_negMult = 1;
	_anchorB = 0;
	_afterBond = -1;
	_threshold = 0.80;
	_increment = 5;
	_useTarget = true;
	_usingWhack = false;
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
	std::cout << "| 3. Performing SVD... " << std::flush;
	Timer timer;
	
	if (_svd)
	{
		delete _svd;
		_svd = NULL;
	}

	_svd = new SVDBond(_bondEffects, _bonds, _atoms);
	_svd->performSVD(&_bbCCs);
	
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

		reflex();
		createAtomTargets();
		scanBondParams();
		
		createClustering();

		/* Testing the SVD stuff */
		svd();
		
		std::vector<ParamBandPtr> extras;

		bool reduceShift = false;
		bool success = false;
		
		for (int j = 0; j < 1; j++)
		{
			std::random_shuffle(_paramBands.begin(), _paramBands.end());

			std::cout << "| " << j + 4 << ". Refining bond clusters... " 
			<< std::flush;
			Timer timer;
			int limit = 5;

			NelderMeadPtr nelder = NelderMeadPtr(new RefinementNelderMead());
			nelder->setCycles(200);
			nelder->setVerbose(true);	

			nelder->setEvaluationFunction(getScore, this);

			_svd->addToStrategy(nelder);

			nelder->refine();

			double val = (1 + getScore(this)) * 100.;
			val = nelder->improvement();

			if (val < 0.5 && !_direct)
			{
				reduceShift = true;
			}

			if (nelder->didChange())
			{
				std::cout << std::setw(3) << val << 
				"% improved. ... done. ";
				
				if (val > 0.5 || _direct)
				{
					timer.quickReport();
					std::cout << std::endl;
					success = true;
					break;
				}
			}
			else
			{
				std::cout << " not improved.   ... done. ";
			}

			timer.quickReport();
			std::cout << std::endl;
			
			if (_paramBands.size() < limit)
			{
				break;
			}
			
		}
		
		std::cout << "---------------------------------------------------------"
		<< std::endl;

		if (reduceShift)
		{
			_shift *= 0.9;
		}

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
//	_paramBands.clear();
	_bbCCs.clear();
	_prepared = false;
}

AtomTarget FlexLocal::currentAtomValues()
{
	AtomTarget targ;

	for (int i = 0; i < _atoms.size(); i++)
	{
		AtomPtr a = _atoms[i];
		double bf = a->getBFactor();

		targ[a] = bf;
	}

	return targ;
}

void FlexLocal::createAtomTargets()
{
	int res = _polymer->getAnchor();
	AtomPtr ca = _polymer->getMonomer(res)->findAtom("CA");
	
	AnchorPtr anchor = _polymer->getAnchorModel();
	_anchorB = anchor->getBFactor();
	
	std::vector<double> all;
	double count = 0;
	_startB = 0;
	
	/* First, choose which atoms to worry about, and find their
	 * starting B factors. */
	for (int i = 0; i < _polymer->atomCount(); i++)
	{
		AtomPtr a = _polymer->atom(i);

		ModelPtr m = a->getModel();
		
		if (!m->isBond())
		{
			continue;
		}
		
		BondPtr b = ToBondPtr(a->getModel());

		if (_usingWhack && !b->getWhack())
		{
			continue;
		}
		
		if (!_usingWhack && 
		   ((!b->getRefineFlexibility() || !b->isNotJustForHydrogens()
			|| !b->isUsingTorsion())))
		{
			continue;
		}
		
		if (a->getResidueNum() == _polymer->getAnchor() && _afterBond < 0)
		{
			_afterBond = _bonds.size();
		}
		
		if (a->getAtomName() != "CA")
		{
			continue;
		}

		_bonds.push_back(b);

		double ibf = a->getInitialBFactor() - ca->getInitialBFactor();
		ibf -= (a->getBFactor() - ca->getBFactor());

		if (_useTarget)
		{
			ibf = a->getTargetB();
		}

		double bf = a->getBFactor();

		_atoms.push_back(a);
		_atomTargets[a] = ibf;
		_atomOriginal[a] = bf;
		
		_startB += bf;
		count++;

		all.push_back(ibf);	
	}
	
	_startB /= count;
	
	if (_useTarget)
	{
		double mean_b = mean(all);
		double stdev_b = standard_deviation(all);

		for (int i = 0; i < _atoms.size(); i++)
		{
			AtomPtr a = _atoms[i];
			_atomTargets[a] -= mean_b;
			_atomTargets[a] += stdev_b / 2;
		}
	}
}	

double FlexLocal::bondAtomCorrel(BondPtr b)
{
	std::vector<double> xs, ys;

	for (int i = 0; i < _atoms.size(); i++)
	{
		AtomPtr a = _atoms[i];
		
		double flex = _bondEffects[b][a];
		double target = targetForAtom(a);
		xs.push_back(flex);
		ys.push_back(target);
	}

	double correl = correlation(xs, ys);
	
	return correl;
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

double FlexLocal::directSimilarity()
{
	std::vector<double> xs, ys;

	for (int i = 0; i < _atoms.size(); i++)
	{
		double target = targetForAtom(_atoms[i]);
		double actual = actualAtomChange(_atoms[i]);
		
		xs.push_back(target);
		ys.push_back(actual);
	}
	
	double rfactor = r_factor(xs, ys);
	
	return rfactor;
}

double FlexLocal::sgetTotalBChange(void *object)
{
	FlexLocal *local = static_cast<FlexLocal *>(object);
	return local->getTotalBChange();
}

double FlexLocal::getTotalB()
{
	double sum = 0;
	double count = 0;

	for (int i = 0; i < _atoms.size(); i++)
	{
		AtomPtr a = _atoms[i];
		double bf = a->getBFactor();

		sum += bf;
		count++;
	}

	sum /= count;

	return sum;
}

double FlexLocal::getTotalBChange()
{
	double sum = getTotalB();
	double diff = _startB - sum;
	
	return fabs(diff);
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
	
	if (local->_direct)
	{
		double score = local->_flexGlobal->score(local->_flexGlobal);
		return score;
	}

//	double score = local->directSimilarity();
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
	
	AtomPtr check = _polymer->findAtoms("CA", 70)[0];
	std::map<int, double> sumsPerAtom;

	for (int i = 0; i < _bonds.size(); i++)
	{
		BondPtr b = _bonds[i];
		WhackPtr w = _bonds[i]->getWhack();
		
		if (!w && _usingWhack)
		{
			continue;
		}
		
		double k = getBondParam(b);

		for (int r = 0; r < 1; r++)
		{
			double add = _shift * _negMult;
			double num = (r == 0) ? i : i + 0.5;

//			std::cout << "Scanning " << b->shortDesc() << std::endl;
			setBondParam(b, k + add);
			
			for (int j = 0; j < _atoms.size(); j++)
			{
				double nAtom = _atoms[j]->getResidueNum();
				double bnow = _atoms[j]->getBFactor();
				double bthen = _atomOriginal[_atoms[j]];

				double diff = bnow - bthen;
				
				if (r == 0)
				{
					if (_bondEffects.count(b) == 0)
					{
						_bondEffects[b] = AtomTarget();
					}

					_bondEffects[b][_atoms[j]] = diff;
				}
				
				int intAtom = _atoms[j]->getResidueNum();

				if (sumsPerAtom.count(intAtom) == 0)
				{
					sumsPerAtom[intAtom] = 0;
				}
				
				sumsPerAtom[intAtom] += diff;

				csv->addEntry(3, (double)nAtom, (double)count, diff);

			}

			setBondParam(b, k);
		}

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

void FlexLocal::setWhacking(bool whack)
{
	_usingWhack = true;
	_getter = Whack::getWhack;
	_setter = Whack::setWhack;
}

double FlexLocal::getBondParam(BondPtr b)
{
	double k = 0;

	if (_usingWhack)
	{
		WhackPtr w = b->getWhack();
		k = (*_getter)(&*w);
	}
	else
	{
		k = (*_getter)(&*b);
	}
	
	return k;
}

void FlexLocal::setBondParam(BondPtr b, double k)
{
	if (_usingWhack)
	{
		WhackPtr w = b->getWhack();
		(*_setter)(&*w, k);
		
		if (*_setter == Whack::setWhack)
		{
			propagateWhack();
		}
	}
	else
	{
		(*_setter)(&*b, k);
		b->propagateChange();
	}
}

