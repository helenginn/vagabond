//
//  BoneDensity.cpp
//  vagabond
//
//  Created by Helen Ginn on 31/03/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "BoneDensity.h"
#include "maths.h"
#include "AtomGroup.h"
#include "Crystal.h"
#include "Monomer.h"
#include "Polymer.h"
#include "Shouter.h"
#include "CSV.h"
#include "polyfit.h"

BoneDensity::BoneDensity()
{
	_start = -1;
	_end = -1;
}

bool BoneDensity::liesInRange(int i)
{
	if (_start < 0 || _end < 0)
	{
		return true;
	}
	
	else return (i >= _start && i <= _end);
}

void BoneDensity::createRefinementStrategies()
{
	int anchor = _polymer->getAnchor();
	int stage = 0;
	
	while (stage < 2)
	{
		int end = (stage == 0 ? 0 : _polymer->monomerCount());
		int skip = (stage == 0) ? -1 : 1;
		int start = (stage == 0) ? anchor - 2 : anchor + 1;
		
		int refineStart = start;
		RefinementType lastR = RefinementFine;

		for (int i = start; i != end; i += skip)
		{
			if (!liesInRange(i))	
			{
				continue;
			}

			if (_summaryMap.count(i) == 0)
			{
				continue;	
			}
			
			double value = _summaryMap[i];
			if (value < 0 && abs(i - anchor) < 4)
			{
				std::cout << "Ignoring resi " << i << " - too close to "\
				"anchor." << std::endl;
				continue;
			}

			RefinementType rType = ((value > 0) ? RefinementFine :
			                        RefinementModelRMSDZero);
			lastR = ((value < 0) ? RefinementFine :
			         RefinementModelRMSDZero);

			int refineEnd = i + skip * 3;
			
			BackboneInstruction inst;
			inst.startRes = refineStart;
			inst.endRes = refineEnd;
			inst.rType = rType;
			_instructions.push_back(inst);
			
			refineStart = i - skip * 3;
		}
		
		BackboneInstruction inst;
		inst.startRes = refineStart;
		inst.endRes = end;
		inst.rType = lastR;
		_instructions.push_back(inst);

		stage++;
	}
}

void BoneDensity::analyse()
{
	/* Check everything is set that needs setting */
	validate();
	_densityMap.clear();
	_summaryMap.clear();
	_instructions.clear();

	std::cout << "*************************" << std::endl;
	std::cout << "Backbone density analysis" << std::endl;
	perMonomerScores();
	
	std::cout << "Finding inflection points..." << std::endl;
	findInflections();
	
	std::cout << "Found " << _summaryMap.size() << " inflection points." << std::endl;
	
	createRefinementStrategies();

}

BackboneState BoneDensity::stateOfBackbone(int start, int end)
{
	if (end < start)
	{
		int tmp = end;
		end = start;
		start = tmp;
	}

	int anchor = _polymer->getAnchor();
	std::vector<double> xs, ys;
	
	/* This shouldn't really ever be called */
	if (anchor > start && anchor < end)
	{
		return BackboneStraddlesAnchor;
	}

	/* Prepare data for polyfit */
	for (int i = start; i <= end; i++)
	{
		if ((_densityMap.count(i) == 0))
		{
			continue;	
		}

		xs.push_back(i);
		ys.push_back(_densityMap[i]);
	}
	
	if (xs.size() < 2)
	{
		std::cout << "Not enough data points" << std::endl;
		return BackboneUnassigned;
	}

	/* Fit a quadratic polynomial */
	std::vector<double> poly = polyfit(xs, ys, 2);

	double inflect = -poly[1] / (2 * poly[2]);
	
	std::cout << "Inflection point at " << inflect << std::endl;

	bool inRange = (inflect >= start && inflect <= end);
	
	/* This is a dubious situation so we make no verdict */
	if (inRange) 
	{
		std::cout << "Inflection point in range" << std::endl;
		return BackboneUnassigned;
	}
	
	std::vector<double> linear = polyfit(xs, ys, 1);

	bool positive = (linear[1] > 0);
	
	/* If we are on the N-terminal side of the anchor point,
	* 	we expect the reverse situation. */
	if (start < anchor)
	{
		positive = !positive;
	}

	/* Some kind of threshold needs to be reached */
	double stdev = standard_deviation(ys);
	double ave = mean(ys);
	double standerr = stdev / ave;
	
	std::cout << "Standard error: " << standerr << std::endl;
	if (standerr < 0.1)
	{
		std::cout << "Error is not high enough" << std::endl;
		return BackboneUnassigned;
	}
	
	if (!positive) 
	{
		std::cout << "Backbone is contracting" << std::endl;
		return BackboneContracting;
	}

	std::cout << "Backbone is expanding" << std::endl;
	return BackboneExpanding;
}

void BoneDensity::findInflections()
{
	std::vector<double> xs, ys;
	int window = 4;
	int direction = 0;
	int run = 0;
	int prevDir = 0;
	double biggestVal = 0;
	int biggestResidue = 0;
	
	for (int i = window; i < _polymer->monomerCount() - window; i++)
	{
		if (!liesInRange(i))	
		{
			continue;
		}

		for (int j = -window; j <= window; j++)
		{
			if ((_densityMap.count(i + j) == 0))
			{
				continue;	
			}

			xs.push_back(i + j);
			ys.push_back(_densityMap[i + j]);
		}
		
		if (xs.size() < 2)
		{
			continue;	
		}

		std::vector<double> poly = polyfit(xs, ys, 2);
		
		double inflect = -poly[1] / (2 * poly[2]);
		bool inRange = (inflect > i - window && inflect <= i + window);
		bool positive = (poly[2] > 0);
		
		if (inRange)
		{
			if (run == 0)
			{
				/* Inflection point needs to be opposite of the last one */
				if ((direction == 0) || (((direction != 0) && 
				                          ((direction < 0 && positive) || 
				                           (direction > 0 && !positive)))))
				{
					prevDir = positive ? 1 : -1;
					run++;
					
					if (fabs(poly[2]) > biggestVal)
					{
						biggestVal = fabs(poly[2]);
						biggestResidue = i;
					}
				}
			}
			else if ((prevDir < 0 && !positive) || (prevDir > 0 && positive))
			{
				/* Collapse me later with repeat above */
				run++;
				if (fabs(poly[2]) > biggestVal)
				{
					biggestVal = fabs(poly[2]);
					biggestResidue = i;
				}
			}
			
			if (run >= 4)
			{
				direction = prevDir;
				_summaryMap[biggestResidue] = biggestVal * direction;

				run = 0;
				biggestResidue = 0;
				biggestVal = 0;
				
			}
		}
		else
		{
			run = 0;
			
			if (direction == 0)
			{
				prevDir = 0;
			}
		}

		xs.clear();
		ys.clear();
	}
}

void BoneDensity::perMonomerScores()
{	
	CSVPtr csv = CSVPtr(new CSV(2, "resnum", "gradient"));
	double sum = 0;
	double count = 0;
	
	for (int i = 0; i < _polymer->monomerCount(); i++)
	{
		if (!liesInRange(i))	
		{
			continue;
		}

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
		_densityMap[i] = gradient;
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
	}
	
	for (DensityScoreMap::iterator it = _densityMap.begin();
	     it != _densityMap.end(); it++)
	{
		it->second /= sum;
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
	
	if (_densityMap.count(anchor) != 0)
	{
		std::cout << "Anchor (" << anchor << ") for chain " << _polymer->getChainID() << " has score: " << _densityMap[anchor] << std::endl;
	}
	
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

void BoneDensity::setRange(int start, int end)
{
	if (end < start)
	{
		int tmp = end;
		end = start;
		start = tmp;
	}
	
	_start = start;
	_end = end;
}

