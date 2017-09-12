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
				if (rType != RefinementModelRMSD)
				{
					break;
				}

				ScoreType scoreType = ScoreTypeModelRMSDZero;
				for (int l = 0; l < 1; l++)
				{
					setupNelderMead();
					bond->setBlocked(true);
					addMagicAxis(bond, deg2rad(20.0), deg2rad(2.0));

					if (scoreType == ScoreTypeModelRMSD)
					{
						addDampening(bond, 0.1, 0.1);
						addTorsionBlur(bond, 0.1, 0.1);
					}

					setJobName("magic_axis_" + i_to_str(resNum) + "_" + bond->shortDesc());
					addSampledAtoms(shared_from_this());
					setSilent();
					setScoreType(scoreType);
					sample();

					if (scoreType == ScoreTypeModelRMSD)
					{
						scoreType = ScoreTypeModelRMSDZero;
					}
					else
					{
						scoreType = ScoreTypeModelRMSD;

					}

					bond->resetAxis();

				}

				setupNelderMead();
				setupTorsionSet(bond, k, 5, resNum, deg2rad(0.2), deg2rad(0.01));
				addSampledAtoms(shared_from_this());
				setScoreType(ScoreTypeModelPos);
				setJobName("model_pos_" + i_to_str(resNum) + "_" + bond->shortDesc());
				sample();

			}

			if (rType == RefinementModelRMSD)
			{
				continue;
			}

			for (int k = 0; k < groups; k++)
			{
				setupNelderMead();
				setupTorsionSet(bond, k, 5, resNum, 2.5, 0.1);
				setCrystal(target);
				sample();

				setupNelderMead();
				setupTorsionSet(bond, k, 5, resNum, 0.5, 0.1);
				addDampening(bond, 0.2, 0.1);
				addTorsionBlur(bond, 0.3, 0.1);
				setCrystal(target);
				sample();
			}


			bond->setActiveGroup(0);

		}

	}
}
