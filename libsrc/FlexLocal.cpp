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
	_shift = 0.10;
}

void FlexLocal::createAtomTargets()
{
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
		
		_atoms.push_back(a);
		_atomTargets[a] = bf;
		
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
				double bthen = _atomTargets[_atoms[j]];

				double diff = bnow - bthen;
				
				if (fabs(Bond::getKick(&*b)) < 1e-6)
				{
					diff = 0;
				}
				
				double grad = diff;

				csv->addEntry(3, (double)nAtom, (double)num, grad);

			}

			Bond::setKick(&*b, k);
			b->propagateChange();

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
}
