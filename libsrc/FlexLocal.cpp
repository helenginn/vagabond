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
	_shift = 0.05;
	_run = 0;
	_window = 10;
	_direct = 0;
	_negMult = 1;
	_anchorB = 0;
	_afterBond = -1;
	_threshold = 0.80;
	_increment = 8;
	_useTarget = false;
	_usingWhack = false;
	_getter = Bond::getKick;
	_flexGlobal = NULL;
	_setter = Bond::setKick;
}

FlexLocal::~FlexLocal()
{
	if (_flexGlobal)
	{
		delete _flexGlobal;
	}
}

void FlexLocal::toElectronDensity()
{
	_flexGlobal = new FlexGlobal();
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	AtomGroupPtr allBackbone = _polymer->getAllBackbone();
	_flexGlobal->setAtomGroup(allBackbone);
	_flexGlobal->setCrystal(crystal);
	_flexGlobal->matchElectronDensity();
	_direct = 1;
}

void FlexLocal::setPolymer(PolymerPtr pol, double shift)
{
	_polymer = pol;
	_shift = shift;
}

void FlexLocal::refineAnchor()
{
	std::cout << "---------------------------------------------------------"
	<< std::endl;
	std::cout << "|  Refining whole molecule for chain " 
	<< _polymer->getChainID();
	std::cout << std::endl;
	std::cout << "---------------------------------------------------------"
	<< std::endl;

	reflex();
	createAtomTargets();
	
	std::cout << "| 1. Refining libration...   " << std::flush;

	Timer timer;
	AnchorPtr anchor = ToAnchorPtr(_polymer->getAnchorModel());

	NelderMeadPtr nelder = NelderMeadPtr(new RefinementNelderMead());
	nelder->setCycles(24);
	nelder->setSilent();	
	
	nelder->setEvaluationFunction(getScore, this);
//	anchor->addTranslationParameters(nelder, 0.2);
	anchor->addLibrationParameters(nelder, 0.2);
	nelder->refine();

	double val = (1 - getScore(this)) * 100.;

	if (nelder->didChange())
	{
		std::cout << std::setw(8) << val << "% improved. ... done. ";
	}
	else
	{
		std::cout << " not improved.";
	}

	timer.quickReport();
	std::cout << std::endl;
	
	std::cout << "---------------------------------------------------------"
	<< std::endl;
}

void FlexLocal::refine()
{
	for (int i = 0; i < 2; i++)
	{
		scanBondParams();
		createClustering();
		reorganiseBondOrder();
		chooseBestDifferenceThreshold();
		std::vector<ParamBandPtr> extras;

		bool reduceShift = false;

		for (int j = 0; j < 5; j++)
		{
			std::cout << "| " << j + 5 << ". Refining bond clusters... " 
			<< std::flush;
			Timer timer;
			int limit = 5;

			NelderMeadPtr nelder = NelderMeadPtr(new RefinementNelderMead());
			nelder->setCycles(30);
			nelder->setSilent();	

			nelder->setEvaluationFunction(getScore, this);

			for (int i = 0; i < _paramBands.size() && i < limit; i++)
			{
				double dir = 1;
				ParamBandPtr band = _paramBands[i];
				nelder->addParameter(&*band, ParamBand::getGlobalParam,
				                     ParamBand::setGlobalParam, _shift * dir, 
				                     _shift / 20, "w" + i_to_str(i));
				
				if (_usingWhack)
				{
					ParamBandPtr second;
					second = ParamBandPtr(new ParamBand(*band));

					second->setPrivateGetter(Whack::getKick);
					second->setPrivateSetter(Whack::setKick);
					second->prepare();
					extras.push_back(second);

					nelder->addParameter(&*second, ParamBand::getGlobalParam,
					                     ParamBand::setGlobalParam,
					                     _shift * dir, 
					                     _shift / 20, "k" + i_to_str(i));
				}
			}

			nelder->refine();

			double val = (1 + getScore(this)) * 100.;
			val = nelder->improvement();

			if (val < 0.5 && !_direct)
			{
				reduceShift = true;
			}

			if (nelder->didChange())
			{
				std::cout << std::setw(5) << val << 
				"% improved. ... done. ";
				
				if (val > 0.5 || _direct)
				{
					timer.quickReport();
					std::cout << std::endl;
					break;
				}
			}
			else
			{
				std::cout << " not improved.   ... done. ";
			}

			timer.quickReport();
			std::cout << std::endl;
			std::random_shuffle(_paramBands.begin(), _paramBands.end());
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
	_paramBands.clear();
	_bbCCs.clear();
}

int getHIndex(std::map<int, int> pairs)
{
	int h = 0;
	int changed = 1;
	while (changed)
	{
		changed = 0;
		for (int i = 0; i < pairs.size(); i++)
		{
			int members = pairs[i];
			if (members > h)
			{
				h++;
				pairs[i] = 0;
				changed = 1;
			}
		}
	}
	
	return h;
}

std::map<int, int> FlexLocal::getClusterMembership(double threshold)
{
	std::map<int, int> result;
	
	int current = 0;
	result[current] = 0;
	
	for (int i = 0; i < _b2bDiffs.size(); i++)
	{
		double diff = _b2bDiffs[i];
		
		if (diff < threshold)
		{
			result[current]++;
		}
		else
		{
			current++;
			result[current] = 0;
		}
	}
	
	return result;
}

void FlexLocal::chooseBestDifferenceThreshold()
{
	Timer timer;
	std::cout << "| 4. Subdividing into clusters..." << std::flush;
	/* t threshold for a difference between two consecutive re-ordered
	 * bonds to cause a change of cluster */
	
	int best_h = 0;
	int total = 0;
	double best_t = 0;

	for (double t = 0; t < 1.0; t += 0.05)
	{
		std::map<int, int> result;
		result = getClusterMembership(t);
		int h = getHIndex(result);
		
		if (best_h <= h)
		{
			best_h = h;
			best_t = t;
			total = result.size();
		}
	}
	
	std::cout << " " << total << " clusters.  ... done. ";
	
	timer.quickReport();
	std::cout << std::endl;

	CSVPtr csv = CSVPtr(new CSV(2, "id_bondnum", "id_cluster"));
	int current = 0;

	ParamBandPtr band;

	double sum = 0;

	for (int i = 0; i < _b2bDiffs.size(); i++)
	{
		if (_b2bDiffs[i] > best_t || i == 0)
		{
			if (band)
			{
				band->prepare();
				sum /= band->objectCount();
				band->setTag(sum);
				sum = 0;
			}

			band = ParamBandPtr(new ParamBand());
			band->setPrivateGetter(_getter);
			band->setPrivateSetter(_setter);

			_paramBands.push_back(band);
			current++;
		}

		BondPtr bond = _bonds[_reorderedBonds[i]];
		sum += bondAtomCorrel(bond);
		
		if (_usingWhack)
		{
			WhackPtr whack = bond->getWhack();
			
			if (whack)
			{
				band->addObject(&*whack, 1);
			}
			
		}
		else
		{
			band->addObject(&*bond, 1);
		}
		
		_bondClusterIds[bond] = current;
		csv->addEntry(2, (double)(_reorderedBonds[i]), (double)current);
	}
	
	band->prepare();
	
	std::random_shuffle(_paramBands.begin(), _paramBands.end());
	
	csv->setSubDirectory("local_flex");
	csv->writeToFile(_polymer->getChainID() + "bond_cluster_ids.csv");

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
	_anchorB = ca->getTargetB();
	
	if (_anchorB != _anchorB || !_useTarget)
	{
		_anchorB = 0;
	}
	
	std::vector<double> all;
	
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
		
		/* Do one side or the other of the Ramachandran angle */
		if (_negMult > 0 && b->getMinor()->getAtomName() == "CA")
		{
			_bonds.push_back(b);
		}
		else if (_negMult < 0 && b->getMajor()->getAtomName() == "CA")
		{
			_bonds.push_back(b);
		}

		/* We don't want to include this as a target atom unless it's a
		 *  C-alpha atom */
		if (a->getAtomName() != "CA")
		{
			continue;
		}
		
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
		all.push_back(ibf);	
	}
	
	if (_useTarget)
	{
		std::sort(all.begin(), all.end(), std::less<double>());
		int aim = all.size() / 20;
		
		if (aim == 0)
		{
			aim = 1;
		}

		double sum = 0;
		
		for (int i = 0; i < aim && i < all.size(); i++)
		{
			sum += all[i];
		}
		
		sum /= (double)aim;
		_anchorB = sum;
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

void FlexLocal::reorganiseBondOrder()
{
	Timer timer;
	std::cout << "| 3. Reorganising bonds into similar groups... " << std::flush;
	
	std::vector<int> orderedResults;
	
	orderedResults.push_back(_degrees[0].index); 
	_degrees.erase(_degrees.begin());
	int parent = 0;

	while (true)
	{
		int added = 0;

		BondPtr pBond = _bonds[orderedResults[parent]];

		/* Go through rest of the array (as they're already degree-ordered)
		 * and queue the ones with a connection to bond no. current. */
		for (int i = 0; i < _degrees.size(); i++)
		{
			BondPtr posChild = _bonds[_degrees[i].index];

			/* Connection */
			if (_bbCCs[pBond][posChild] > _threshold)
			{
				orderedResults.push_back(_degrees[i].index);
				_degrees.erase(_degrees.begin() + i);
				added++;
				i--;
			}
		}
		
		if (orderedResults.size() == _bonds.size())
		{
			break;
		}

		parent++;
		
		if (parent >= orderedResults.size())
		{
			orderedResults.push_back(_degrees[0].index);
			_degrees.erase(_degrees.begin());
		}
	}
	
	std::cout << " ... done. ";
	timer.quickReport();
	std::cout << std::endl;
	
	CSVPtr csv = CSVPtr(new CSV(3, "pbond_i", "pbond_j", "pbondcc"));

	for (int i = 0; i < orderedResults.size(); i++)
	{
		double bondnum = orderedResults[i];
		BondPtr bi = _bonds[orderedResults[i]];

		for (int j = 0; j < orderedResults.size(); j++)
		{
			BondPtr bj = _bonds[orderedResults[j]];
			
			double cc = _bbCCs[bi][bj];

			csv->addEntry(4, (double)i, (double)j, cc);
		}
	}

	csv->setSubDirectory("local_flex");
	csv->writeToFile(_polymer->getChainID() + "_after_bond_matrix.csv");
	
	for (int i = 0; i < orderedResults.size() - 1; i++)
	{
		double bond_num = orderedResults[i];
		BondPtr ba = _bonds[orderedResults[i]];
		BondPtr bb = _bonds[orderedResults[i + 1]];

		double diffs = 0;
		int count = 0;

		for (int j = i + 1; j < orderedResults.size() && j < i + 30; j++)
		{
			BondPtr bj = _bonds[orderedResults[j]];
			double cc_a = _bbCCs[bj][ba];
			double cc_b = _bbCCs[bj][bb];
			
			double diff = (cc_b - cc_a) * (cc_b - cc_a);
			diffs += diff;
			count++;
		}
		
		diffs /= (double)count;
		
		_b2bDiffs.push_back(diffs);
	}
	
	_reorderedBonds = orderedResults;
}

void FlexLocal::createClustering()
{
	CSVPtr csv = CSVPtr(new CSV(3, "bond_i", "bond_j", "bondcc"));

	int anchor = _polymer->getAnchor();
	double sumCC = 0;
	double count = 0;

	Timer timer;
	std::cout << "| 2. Calculating bond-bond correlation pairs... " << std::flush;
	/* Once all the clusters are made, give them pair-wise correlations */
	for (int i = 0; i < _bonds.size(); i++)
	{
		BondPtr bi = _bonds[i];
		
		int degree = 0;

		for (int j = 0; j < _bonds.size(); j++)
		{
			BondPtr bj = _bonds[j]; // giggle

			double cc = 0;
			if (bi < bj)
			{
				cc = bondRelationship(bi, bj);
			}
			else
			{
				cc = _bbCCs[bi][bj];
			}
			
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

double FlexLocal::getScore(void *object)
{
	FlexLocal *local = static_cast<FlexLocal *>(object);
	local->_polymer->propagateChange();
	
	if (local->_direct)
	{
		double score = local->_flexGlobal->score(local->_flexGlobal);
		return score;
	}

	double score = local->directSimilarity();
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
	reflex.setPieceCount(3);
	reflex.calculate();
	std::cout << "   ... done. ";
	timer.quickReport();
	std::cout << std::endl;
}

void FlexLocal::scanBondParams()
{
	std::cout << "---------------------------------------------------------"
	<< std::endl;
	std::cout << "|  Refining flexibility for chain " 
	<< _polymer->getChainID() << " (cycle " << _run << ")";
	std::cout << std::endl;
	std::cout << "---------------------------------------------------------"
	<< std::endl;
	
	ExplicitModelPtr model = _polymer->getAnchorModel();
	int samples = model->getFinalPositions().size();

	reflex();
	createAtomTargets();

	Options::setNSamples(40);
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
		plotMap["filename"] = "atom_targets_" + _polymer->getGraphName();
		plotMap["height"] = "700";
		plotMap["width"] = "1200";
		plotMap["xHeader0"] = "atom";
		plotMap["yHeader0"] = "target";
		plotMap["colour0"] = "black";
		plotMap["xTitle0"] = "Residue number";
		plotMap["yTitle0"] = "B factor target";
		plotMap["style0"] = "line";

		atomtarg->plotPNG(plotMap);

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
			double add = _shift;
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
	plotMap["filename"] = "bond_scan_" + _polymer->getGraphName();
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

	Options::setNSamples(samples);
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

