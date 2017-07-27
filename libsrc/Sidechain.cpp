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

void Sidechain::refine(CrystalPtr target)
{
	int resNum = getMonomer()->getResidueNum();

	for (int i = 0; i < bondCount(); i++)
	{
		if (bond(i)->isNotJustForHydrogens() && bond(i)->isUsingTorsion())
		{
			addTorsion(bond(i), deg2rad(60), deg2rad(1.0));
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