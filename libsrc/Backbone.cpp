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

			if (!bond->isNotJustForHydrogens() || bond->isFixed()
				|| !bond->isUsingTorsion())
			{
				continue;
			}

			int groups = bond->downstreamAtomGroupCount();

			for (int k = 0; k < groups; k++)
			{
				if (rType != RefinementModelRMSD)
				{
					break;
				}

				ScoreType scoreType = ScoreTypeModelRMSD;
				for (int l = 0; l < 3; l++)
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
					addSampledCAs(getPolymer(), resNum, resNum + 24);
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
				setJobName("model_pos_" + bond->shortDesc());
				addRamachandranAngles(getPolymer(), resNum, resNum + 3);
				addSampledCAs(getPolymer(), resNum, resNum + 3);
				setScoreType(ScoreTypeModelPos);
				sample();

			}

			if (rType == RefinementModelRMSD)
			{
				continue;
			}

			for (int k = 0; k < groups; k++)
			{
				bool shouldBlur = (rType == RefinementFineBlur);

				if (shouldBlur)
				{
					setupNelderMead();
					setupTorsionSet(bond, k, 5, resNum, 0.05, 0.02, true);
					setCrystal(target);
					sample();
				}

				setupNelderMead();
				setupTorsionSet(bond, k, 4, resNum, 0.2, 0.1);

				setCrystal(target);
				sample();
			}


			bond->setActiveGroup(0);

		}

	}
}
