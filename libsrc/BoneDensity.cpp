//
//  BoneDensity.cpp
//  vagabond
//
//  Created by Helen Ginn on 31/03/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "BoneDensity.h"
#include "AtomGroup.h"
#include "Crystal.h"
#include "Monomer.h"
#include "Polymer.h"
#include "Shouter.h"
#include "CSV.h"

typedef std::map<int, double> DensityScoreMap;

BoneDensity::BoneDensity()
{

}


void BoneDensity::analyse()
{
	/* Check everything is set that needs setting */
	validate();

	std::cout << "Backbone density analysis" << std::endl;

	/* Score gradient of every monomer */
	
	DensityScoreMap densityMap;
	CSVPtr csv = CSVPtr(new CSV(2, "resnum", "gradient"));
	double sum = 0;
	double count = 0;
	
	for (int i = 0; i < _polymer->monomerCount(); i++)
	{
		AtomGroupPtr trio = AtomGroupPtr(new AtomGroup());
			
		for (int j = -1; j < 2; j++)
		{
			MonomerPtr monomer = _polymer->getMonomer(i + j);
			if (!monomer)
			{
				continue;	
			}
			
			trio->addAtomsFrom(monomer);
		}
		
		if (trio->atomCount() == 0)
		{
			continue;	
		}
		
		double gradient = -trio->scoreWithMap(ScoreTypeMultiply, _crystal);
		csv->addEntry(2, (double)i, gradient);
		sum += gradient;
		count++;
	}
	
	sum /= count;
	
	for (int i = 0; i < csv->entryCount(); i++)
	{
		double grad = csv->valueForEntry("gradient", i);
		grad /= sum;
		csv->setValueForEntry(i, "gradient", grad);
		densityMap[i] = grad;
	}
	
	std::cout << "Weighted sum of map voxels for monomer trios." << std::endl;
	
	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = "gradients";
	plotMap["xHeader0"] = "resnum";
	plotMap["yHeader0"] = "gradient";
	plotMap["colour0"] = "black";

	plotMap["xTitle0"] = "Residue number";
	plotMap["yTitle0"] = "Gradient";
	plotMap["style0"] = "line";

	csv->writeToFile("gradients.csv");
	csv->plotPNG(plotMap);

	int anchor = _polymer->getAnchor();
	
	if (densityMap.count(anchor) == 0)
	{
		shout_at_helen("Somehow, the density anchor is missing during/n"\
		               "backbone density analysis.");
	}
	
	std::cout << "Anchor for chain " << _polymer->getChainID() << " has score: " << densityMap[anchor] << std::endl;
}



void BoneDensity::validate()
{
	if (!_crystal)
	{
		shout_at_helen("Trying to analyse backbone density\nwithout a crystal.");	
	}

	if (!_polymer)
	{
		shout_at_helen("Trying to analyse backbone density\nwithout a polymer.");	
	}
}
