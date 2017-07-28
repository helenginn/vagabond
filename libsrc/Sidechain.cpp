//
//  Sidechain.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Sidechain.h"
#include "Sampler.h"
#include "Bond.h"
#include "Atom.h"
#include <iostream>
#include "Monomer.h"
#include "FileReader.h"

void Sidechain::refine(CrystalPtr target, RefinementType rType)
{
	if (rType == RefinementBroad)
	{
		setupGrid();
	}
	else
	{
		setupGrid();
//		setupNelderMead();
	}

	int resNum = getMonomer()->getResidueNum();

	for (int i = 0; i < bondCount(); i++)
	{
		if (bond(i)->isNotJustForHydrogens() && bond(i)->isUsingTorsion())
		{
			if (rType == RefinementBroad)
			{
				addTorsion(bond(i), deg2rad(30), deg2rad(2.0));
			}
			else if (rType == RefinementFine)
			{
				addTorsion(bond(i), deg2rad(20.0), deg2rad(0.2));
				addTorsionBlur(bond(i), deg2rad(3.0), deg2rad(0.2));
			}
		}
	}

	setJobName("sample_torsion_" + i_to_str(resNum));

	if (!sampleSize())
	{
		for (int i = 0; i < atomCount(); i++)
		{
			if (atom(i)->getAtomName() != "CA")
			{
				addSampled(atom(i));
			}
		}
	}

	setCrystal(target);
	sample();

}