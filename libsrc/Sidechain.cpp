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
#include "Absolute.h"

void Sidechain::refine(CrystalPtr target, RefinementType rType)
{
	int resNum = getMonomer()->getResidueNum();

	if (!isTied())
	{
		return;
	}

	_timesRefined++;

	const char *atoms[] = {"CB", "OG", "CG", "SD", "CG1", "CG2", "CD", "CE", "NZ", "OG1"};
	std::vector<std::string> atomStrs(atoms, atoms+8);

	for (int i = 0; i < atomStrs.size(); i++)
	{
		std::string atom = atomStrs[i];
		AtomList myAtoms = findAtoms(atom);

		if (!myAtoms.size())
		{
			continue;
		}

		for (int j = 0; j < myAtoms.size(); j++)
		{
			AtomPtr myAtom = myAtoms[j].lock();
			ModelPtr model = myAtom->getModel();

			BondPtr bond = std::static_pointer_cast<Bond>(model);
			std::string majorAtom = bond->getMajor()->getAtomName();

			if (!bond->isNotJustForHydrogens() || bond->isFixed())
			{
				continue;
			}

			int groups = bond->downstreamAtomGroupCount();

			ModelPtr preModel = bond->getParentModel();

			if (bond->isUsingTorsion() && bond->isNotJustForHydrogens())
			{
				if (rType == RefinementBroad)
				{
					for (int k = 0; k < groups; k++)
					{
						setupDoubleTorsion(bond, k, 0, resNum, 360, 8);
						setCrystal(target);
						sample();

						setupDoubleTorsion(bond, k, 0, resNum, 30, 2);
						setCrystal(target);
						sample();
					}
				}
				else
				{
					for (int k = 0; k < groups; k++)
					{
						bond->setActiveGroup(k);

						setupNelderMead();
						reportInDegrees();

						addTorsion(bond, deg2rad(0.2), deg2rad(0.5));

						for (int j = 0; j < bond->downstreamAtomCount(k); j++)
						{
							addSampled(bond->downstreamAtom(k, j));
						}

						for (int j = 0; j < bond->extraTorsionSampleCount(k); j++)
						{
							addSampled(bond->extraTorsionSample(k, j));
						}

						setJobName("torsion_" + majorAtom + "_" + atom + "_g" +
								   i_to_str(k) + "_" + i_to_str(resNum));
						setCrystal(target);
						sample();
					}
					
					bond->setActiveGroup(0);
				}

				if (false && bond->isUsingTorsion() && (rType == RefinementFine))
				{
					for (int k = 0; k < groups; k++)
					{
						bond->setActiveGroup(k);
						setupNelderMead();

						addTorsionNextBlur(bond, 0.2, 0.5);
					//	addTorsionBlur(bond, deg2rad(12.0), deg2rad(0.5));

						for (int j = 0; j < bond->downstreamAtomCount(k); j++)
						{
							addSampled(bond->downstreamAtom(k, j));
						}

						setJobName("blur_" + majorAtom + "_"
								   + atom + "_" + i_to_str(resNum));
						setCrystal(target);
						sample();
					}

					bond->setActiveGroup(0);
				}
			}
		}

	}
}
