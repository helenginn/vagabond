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
#include "RefinementGridSearch.h"
#include "RefinementNelderMead.h"
#include "RefinementLBFGS.h"
#include "ParamBand.h"
#include <map>
#include <iomanip>

FlexLocal::FlexLocal()
{
	_shift = 0.05;
	_run = 0;
	_window = 10;
	_direct = 1;
	_anchorB = 0;
	_afterBond = -1;
	_threshold = 0.80;
	_increment = 6;
}

void FlexLocal::setPolymer(PolymerPtr pol, double shift)
{
	_polymer = pol;
	_shift = shift;
}


void FlexLocal::refine()
{
	scanBondParams();
	createClustering();
	reorganiseBondOrder();
	chooseBestDifferenceThreshold();

	std::cout << "5. Refining parameters of bond clusters..." << std::flush;
	RefinementGridSearchPtr grid;
	grid = RefinementGridSearchPtr(new RefinementGridSearch());
	grid->setCycles(24);
	grid->setSilent();
	grid->setEvaluationFunction(getScore, this);
	int limit = 5;

	for (int i = 0; i < _paramBands.size() && i < limit; i++)
	{
		grid->addParameter(&*_paramBands[i], ParamBand::getGlobalParam,
		                   ParamBand::setGlobalParam, _shift * 2, 
		                   _shift, "k" + i_to_str(i));
	}

	grid->refine();
	
	if (!grid->didChange())
	{
		_shift *= 0.8;
		std::cout << " not improved.\nReducing shift for next time ... done." << std::endl;
	}
	
	_direct = 1;

	NelderMeadPtr nelder;
	nelder = NelderMeadPtr(new RefinementNelderMead());
	nelder->setCycles(24);
	nelder->setSilent();
	nelder->setEvaluationFunction(getScore, this);
	
	for (int i = 0; i < grid->parameterCount(); i++)
	{
		Parameter param = grid->getParamObject(i);
		ParamBand *band = static_cast<ParamBand *>(param.object);
		
		bool changed = grid->didChange(i);
		
		if (!changed)
		{
			continue;
		}
		
		param.other_value /= 20;
		
		nelder->addParameter(param);
	}
	
	double val = (1 - getScore(this)) * 100.;

	if (nelder->parameterCount() > 0)
	{
		nelder->refine();
	}
	
	if (grid->didChange())
	{
		std::cout << val << "% improved. ... done." << std::endl;
	}
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
	std::cout << "4. Subdividing into clusters..." << std::flush;
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
	
	std::cout << " " << total << " clusters.  ... done." << std::endl;
	
	CSVPtr csv = CSVPtr(new CSV(2, "id_bondnum", "id_cluster"));
	int current = 0;

	ParamBandPtr band;

	for (int i = 0; i < _b2bDiffs.size(); i++)
	{
		if (_b2bDiffs[i] > best_t || i == 0)
		{
			if (band)
			{
				band->prepare();
			}

			band = ParamBandPtr(new ParamBand());
			band->setPrivateGetter(Bond::getKick);
			band->setPrivateSetter(Bond::setKick);

			_paramBands.push_back(band);
			current++;
		}

		BondPtr bond = _bonds[_reorderedBonds[i]];
		band->addObject(&*bond, 1);
		
		_bondClusterIds[bond] = current;
		csv->addEntry(2, (double)(_reorderedBonds[i]), (double)current);
	}
	
	band->prepare();
	
	std::random_shuffle(_paramBands.begin(), _paramBands.end());
	
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
	ExplicitModelPtr model = _polymer->getAnchorModel();
	_anchorB = model->getAtom()->getInitialBFactor() - model->getBFactor();
	
	/* First, choose which atoms to worry about, and find their
	 * starting B factors. */
	for (int i = 0; i < _polymer->atomCount(); i++)
	{
		AtomPtr a = _polymer->atom(i);
		double ibf = a->getInitialBFactor();
		double bf = a->getBFactor();
		
		if (!(a->isBackbone() || a->isBackboneAndSidechain()))
		{
			continue;
		}
		
		if (a->getElectronCount() == 1)
		{
			continue;
		}
		
		_atoms.push_back(a);
		_atomTargets[a] = ibf;
		_atomOriginal[a] = bf;
		
		ModelPtr m = a->getModel();
		
		if (!m->isBond())
		{
			continue;
		}
		
		BondPtr b = ToBondPtr(a->getModel());
		
		if (!b->getRefineFlexibility() || !b->isNotJustForHydrogens()
			|| !b->isUsingTorsion())
		{
			continue;
		}
		
		if (a->getResidueNum() == _polymer->getAnchor() && _afterBond < 0)
		{
			_afterBond = _bonds.size();
		}

		_bonds.push_back(b);
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

bool less_than(BondDegree &b1, BondDegree &b2)
{
	return (b2.degree > b1.degree);
}

void FlexLocal::reorganiseBondOrder()
{
	std::cout << "3. Reorganising bonds into similar groups... " << std::flush;
	
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
	
	std::cout << " ... done." << std::endl;
	
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

	/* Degree of each row in a sparse matrix */
	
	int anchor = _polymer->getAnchor();
	double sumCC = 0;
	double count = 0;

	std::cout << "2. Calculating bond-bond correlation pairs... " << std::flush;
	/* Once all the clusters are made, give them pair-wise correlations */
	for (int i = 0; i < _bonds.size(); i++)
	{
		BondPtr bi = _bonds[i];
		
		int degree = 0;

		for (int j = 0; j < _bonds.size(); j++)
		{
			BondPtr bj = _bonds[j]; // giggle

			double cc = 0;
			if (bi->getAtom()->getResidueNum() < anchor &&
			    bj->getAtom()->getResidueNum() > anchor)
			{
				cc = 0;
			}
			else if (bi < bj)
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
	
	std::cout << "... done." << std::endl;

	csv->writeToFile(_polymer->getChainID() + "bond_matrix.csv");
	
	std::sort(_degrees.begin(), _degrees.end(), less_than);

}

double FlexLocal::directSimilarity()
{
	std::vector<double> xs, ys;
	for (int i = 0; i < _atoms.size(); i++)
	{
		double target = _atomTargets[_atoms[i]];
		double original = _atomOriginal[_atoms[i]];
		double current = _atoms[i]->getBFactor();
		double transformed = (target - _anchorB - original) / _increment;
		xs.push_back(transformed);
		ys.push_back(current - original);
	}
	
	double correl = r_factor(xs, ys);
	return correl;
}

double FlexLocal::getScore(void *object)
{
	FlexLocal *local = static_cast<FlexLocal *>(object);
	local->_polymer->propagateChange();
	
	double score = local->directSimilarity();
	return score;
}

void FlexLocal::scanBondParams()
{
	std::cout << "1. Determining atom-flexibility effects... " << std::flush;
	createAtomTargets();
	
	CSVPtr atomtarg = CSVPtr(new CSV(2, "atom", "target"));
	
	for (int i = 0; i < _atoms.size(); i++)
	{
		double num = i / 4 + _polymer->beginMonomer()->first;
		atomtarg->addEntry(2, num, _atomTargets[_atoms[i]]);
	}

	atomtarg->writeToFile("atom_targets_" + _polymer->getChainID() + ".csv");
	
	CSVPtr csv = CSVPtr(new CSV(3, "atom", "kick", "response"));

	for (int i = 0; i < _bonds.size(); i++)
	{
		BondPtr b = _bonds[i];
		double k = Bond::getKick(&*b);

		for (int r = 0; r < 1; r++)
		{
			double add = (r == 0) ? _shift : -_shift;
			double num = (r == 0) ? i : i + 0.5;

			Bond::setKick(&*b, k + add);
			b->propagateChange();

			for (int j = 0; j < _atoms.size(); j++)
			{
				double nAtom = (double)j / 4. +
				(double) _polymer->beginMonomer()->first;
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
				
				double grad = diff;

				csv->addEntry(3, (double)nAtom, (double)num, grad);

			}

			Bond::setKick(&*b, k);
			b->propagateChange();
		}

	}
	
	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = "chain_" + _polymer->getChainID() + "_bondscan";
	plotMap["height"] = "1000";
	plotMap["width"] = "1000";
	plotMap["xHeader0"] = "kick";
	plotMap["yHeader0"] = "atom";
	plotMap["zHeader0"] = "response";

	plotMap["xTitle0"] = "kick";
	plotMap["yTitle0"] = "atom";

	plotMap["style0"] = "heatmap";
	plotMap["stride"] = i_to_str(_bonds.size());
	
	csv->writeToFile(plotMap["filename"] + ".csv");
	csv->plotPNG(plotMap);
	
	std::cout << "   ... done." << std::endl;
}


