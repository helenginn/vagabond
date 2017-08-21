//
//  Backbone.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Backbone.h"
#include "Monomer.h"
#include "Bond.h"
#include "Atom.h"
#include "Absolute.h"
#include "FileReader.h"

void Backbone::refine(CrystalPtr target, RefinementType rType)
{
	if (rType == RefinementBroad)
	{
	//	return;
	}


	if (!isTied())
	{
		return;
	}

	int resNum = getMonomer()->getResidueNum();

	const char *atoms[] = {"CA", "N", "C"};
	std::vector<std::string> atomStrs(atoms, atoms+2);

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

			if (model->getClassName() == "Absolute")
			{
				continue;
			}

			BondPtr bond = std::static_pointer_cast<Bond>(model);
			std::string majorAtom = bond->getMajor()->getAtomName();

			if (!bond->isNotJustForHydrogens() || bond->isFixed())
			{
				continue;
			}

			int groups = bond->downstreamAtomGroupCount();

			if (false && bond->isUsingTorsion() && (rType == RefinementFine))
			{
				for (int k = 0; k < groups; k++)
				{
					bond->setActiveGroup(k);
					setupNelderMead();

					addDampening(bond, 0.2, 0.5);

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

			if (bond->isUsingTorsion() && bond->isNotJustForHydrogens())
			{
				if (rType == RefinementBroad)
				{
					for (int k = 0; k < groups; k++)
					{
						setupGrid();
						setupDoubleTorsion(bond, k, 1, resNum, 360, 8);
						setCrystal(target);
						sample();

						setupGrid();
						setupDoubleTorsion(bond, k, 1, resNum, 16, 1);
						setCrystal(target);
						sample();
					}

				}
				else
				{
					for (int k = 0; k < groups; k++)
					{
						setupNelderMead();
						setupDoubleTorsion(bond, k, 5, resNum, 0.2, 0.1);
						setCrystal(target);
						sample();
					}
					
					bond->setActiveGroup(0);
					
				}
			}
		}
		
	}
}
