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

			ModelPtr preModel = bond->getParentModel();

			int groups = bond->downstreamAtomGroupCount();

			for (int k = 0; k < groups; k++)
			{
				if (rType != RefinementModelOnly)
				{
					break;
				}

				setupNelderMead();
				bond->setBlocked(true);
				addMagicAxis(bond, deg2rad(20.0), deg2rad(1.0));
				setJobName("magic_axis_" + bond->shortDesc());
				addSampledAtoms(shared_from_this());
				setScoreType(ScoreTypeModelRMSD);
				bond->resetAxis();
				sample();
			}

			if (rType == RefinementModelOnly)
			{
				continue;
			}

			for (int k = 0; k < groups; k++)
			{
				setupNelderMead();
				setupTorsionSet(bond, k, 5, resNum, 0.5, 0.1);
		//		setScoreType(ScoreTypeCorrel);

		//		if (rType == RefinementFineBlur)
				{
					addDampening(bond, 0.2, 0.1);
					addTorsionBlur(bond, 0.3, 0.1);
				}
				setCrystal(target);

				std::cout << "  ";
				sample();
			}


			bond->setActiveGroup(0);

		}

	}
}
