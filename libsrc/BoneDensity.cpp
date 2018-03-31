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

BoneDensity::BoneDensity()
{

}


void BoneDensity::analyse()
{
	/* Check everything is set that needs setting */
	validate();
	
	/* Scale entire map to the backbone density */
	AtomGroupPtr allBones = _polymer->getAllBackbone();
	allBones->scoreWithMap(ScoreTypeScaleOnly, _crystal);

	/* Score gradient of every monomer */
	
	CSVPtr csv = CSVPtr(new CSV(2, "resnum", "gradient"));
	
	for (int i = 0; i < _polymer->monomerCount(); i++)
	{
		MonomerPtr monomer = _polymer->getMonomer(i);
		if (!monomer)
		{
			continue;	
		}
		
		double gradient = monomer->scoreWithMap(ScoreTypeMultiply, _crystal);
		csv->addEntry(2, (double)i, gradient);
	}
	
	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = "gradients";
	plotMap["xHeader0"] = "resnum";
	plotMap["yHeader0"] = "gradient";
	plotMap["colour0"] = "black";

	plotMap["xTitle0"] = "Residue number";
	plotMap["yTitle0"] = "Gradient";
	plotMap["style0"] = "line";

	csv->writeToFile("gradients.csv");
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
