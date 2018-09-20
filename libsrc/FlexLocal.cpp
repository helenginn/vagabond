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
#include "Clusterable.h"
#include "Monomer.h"
#include "Anchor.h"
#include "RefinementGridSearch.h"
#include "RefinementNelderMead.h"
#include <map>
#include <iomanip>

FlexLocal::FlexLocal()
{
	_shift = 0.05;
	_run = 0;
	_window = 10;
	_direct = 0;
	_anchorB = 0;
}

void FlexLocal::refine()
{
	scanBondParams();
	createClustering();
	return;

	_direct = 1;

	RefinementGridSearchPtr grid;
	grid = RefinementGridSearchPtr(new RefinementGridSearch());
	grid->setCycles(24);
	grid->setVerbose(true);
	grid->setEvaluationFunction(getScore, this);
	size_t period = (_bonds.size() / 5);

	for (int i = period / 2;
	     i < _bonds.size(); i += period)
	{
		int rnd = random() % period - 5;
		int n = i + rnd;
		
		if (n < 0 || n >= _bonds.size())
		{
			continue;
		}
		
		grid->addParameter(&*_bonds[n], Bond::getKick, Bond::setKick,
		                   _shift * 2, _shift, "k" + i_to_str(n));
	}

	grid->refine();
	
	_direct = 1;

	NelderMeadPtr nelder;
	nelder = NelderMeadPtr(new RefinementNelderMead());
	nelder->setCycles(24);
	nelder->setVerbose(true);
	nelder->setEvaluationFunction(getScore, this);
	
	for (int i = 0; i < grid->parameterCount(); i++)
	{
		bool changed = grid->didChange(i);
		
		if (!changed)
		{
			continue;
		}
		
		Parameter param = grid->getParamObject(i);
		param.other_value /= 20;
		nelder->addParameter(param);
	}
	
	nelder->refine();
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
	_anchorB = model->getBFactor();
	
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

	double correl = correlation(xs, ys);
	
	if (correl < 0) correl = 0;
	
	return correl;
}

double FlexLocal::clusterScore(void *object)
{
	FlexLocal *me = static_cast<FlexLocal *>(object);
	return me->clusterScore();
}

double FlexLocal::clusterScore()
{
	double fx = 0;

	for (int i = 0; i < _bonds.size(); i++)
	{
		BondPtr b = _bonds[i];
		double sum = _clusters[b]->sumContributionToEval();
		fx += sum;
	}
	
	fx /= (double)(_bonds.size() * _bonds.size());
	
	return fx;
}

void FlexLocal::createClustering()
{
	for (int i = 0; i < _bonds.size(); i++)
	{
		BondPtr bi = _bonds[i];
		ClusterablePtr cluster = ClusterablePtr(new Clusterable(_nClusters));
		_clusters[bi] = cluster;
	}
	
	int anchor = _polymer->getAnchor();

	/* Once all the clusters are made, give them pair-wise correlations */
	for (int i = 0; i < _bonds.size() - 1; i++)
	{
		BondPtr bi = _bonds[i];
		ClusterablePtr ci = _clusters[bi];

		for (int j = i + 1; j < _bonds.size(); j++)
		{
			BondPtr bj = _bonds[j]; // giggle

			/* These shouldn't even talk if they span the anchor point
			 * as they are not interconnected. */
			if (bi->getAtom()->getResidueNum() < anchor &&
			    bj->getAtom()->getResidueNum() > anchor)
			{
				continue;
			}

			ClusterablePtr cj = _clusters[bj];

			double cc = bondRelationship(bi, bj);
			
			if (cc != cc)
			{
				continue;
			}
			
			ci->addRelationship(cj, cc);
		}
	}
	
	std::cout << "Starting cluster score: " << clusterScore() << std::endl;
}

double FlexLocal::directSimilarity()
{
	std::vector<double> xs, ys;
	for (int i = 0; i < _atoms.size(); i++)
	{
		double target = _atomTargets[_atoms[i]];
		double original = _atomOriginal[_atoms[i]];
		double current = _atoms[i]->getBFactor();
		double transformed = (target - _anchorB - original) / 5;
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

		for (int r = 0; r < 2; r++)
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
	
	std::cout << "Determined kick-bond effect on atom flexibility."
	<< std::endl;
}


