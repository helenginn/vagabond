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
#include "Monomer.h"
#include "Anchor.h"
#include <map>
#include <iomanip>

FlexLocal::FlexLocal()
{

}

void FlexLocal::scanBondParams()
{
	MonomerPtr anchoredRes = _polymer->getMonomer(_polymer->getAnchor());
	AtomPtr close = anchoredRes->findAtom("CA");
	double base_b = close->getInitialBFactor();
	base_b -= close->getBFactor();

	
	std::map<AtomPtr, double> bFacs;
	std::map<AtomPtr, double> targets;
	std::vector<BondPtr> bonds;
	std::vector<AtomPtr> atoms;

	for (int i = 0; i < _polymer->atomCount(); i++)
	{
		AtomPtr a = _polymer->atom(i);
		double bf = a->getBFactor();
		
		if (!(a->isBackbone() || a->isBackboneAndSidechain()))
		{
			continue;
		}
		
		if (a->getElectronCount() == 1)
		{
			continue;
		}
		
		atoms.push_back(a);
		bFacs[a] = bf;

		double init_b = a->getInitialBFactor();		

		targets[a] = init_b - bf - base_b;
		
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

		bonds.push_back(b);
	}
	
	CSVPtr atomtarg = CSVPtr(new CSV(2, "atom", "target"));
	
	for (int i = 0; i < atoms.size(); i++)
	{
		double num = i / 4 + _polymer->beginMonomer()->first;
		atomtarg->addEntry(2, num, targets[atoms[i]]);
	}

	atomtarg->writeToFile("atom_targets_" + _polymer->getChainID() + ".csv");
	
	double shift = 0.10;
	CSVPtr csv = CSVPtr(new CSV(3, "atom", "kick", "response"));
	CSVPtr bondcc = CSVPtr(new CSV(2, "bond", "cc"));

	for (int i = 0; i < bonds.size(); i++)
	{
		BondPtr b = bonds[i];
		double k = Bond::getKick(&*b);

		for (int r = 0; r < 2; r++)
		{
			std::vector<double> tx, ty;	
			double add = (r == 0) ? shift : -shift;
			double num = (r == 0) ? i : i + 0.5;

			Bond::setKick(&*b, k + add);
			b->propagateChange();

			for (int j = 0; j < atoms.size(); j++)
			{
				double nAtom = (double)j / 4. +
				(double) _polymer->beginMonomer()->first;
				double bnow = atoms[j]->getBFactor();
				double bthen = bFacs[atoms[j]];

				double diff = bnow - bthen;
				
				if (fabs(Bond::getKick(&*b)) < 1e-6)
				{
					diff = 0;
				}
				
				double grad = diff;

				csv->addEntry(3, (double)nAtom, (double)num, grad);

				if (grad > 1e-6)
				{
					tx.push_back(targets[atoms[j]]);
					ty.push_back(grad);
				}

			}

			Bond::setKick(&*b, k);
			b->propagateChange();
			double correl = correlation(tx, ty);

			if (correl != correl)
			{
				correl = 0;
			}

			if (r == 0)
			{
				std::cout << "Done bond " << i << ", " << b->shortDesc();
				
				if (fabs(k) > 1e-6)
				{
					std::cout << " curr. kick = " << std::setprecision(2) << k << std::endl;
				}
				else
				{
					std::cout << std::endl;
				}
			}

			bondcc->addEntry(2, (double)i, correl);
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
	plotMap["stride"] = i_to_str(bonds.size());
	
	csv->writeToFile(plotMap["filename"] + ".csv");
	bondcc->writeToFile(plotMap["filename"] + "_cc_bond.csv");
	csv->plotPNG(plotMap);
}
