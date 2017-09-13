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

	const char *atoms[] = {"N", "CA", "C"};
	std::vector<std::string> atomStrs(atoms, atoms+3);

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

			BondPtr bond = ToBondPtr(model);
			std::string majorAtom = bond->getMajor()->getAtomName();

			if (!bond->isNotJustForHydrogens() || bond->isFixed()
				|| !bond->isUsingTorsion())
			{
				continue;
			}

			int groups = bond->downstreamAtomGroupCount();

			for (int k = 0; k < groups; k++)
			{
				BondPtr nextBond;
				if (bond->getParentModel()->getClassName() == "Bond")
				{
			//		nextBond = ToBondPtr(bond->getParentModel());
				}

				if (rType != RefinementModelRMSD)
				{
					break;
				}

				ScoreType scoreType = ScoreTypeModelRMSDZero;
				for (int l = 0; l < 1; l++)
				{
					if (scoreType == ScoreTypeModelRMSDZero)
					{
						setupGrid();
						addMagicAxisBroad(bond);
						setSilent();
						setJobName("broad_axis_" + i_to_str(resNum) + "_" + bond->shortDesc());
						addSampledBackbone(getPolymer(), resNum, resNum + 20);
						setScoreType(scoreType);
						sample();
					}

					setupNelderMead();
					addMagicAxis(bond, deg2rad(20.0), deg2rad(2.0));

					if (scoreType == ScoreTypeModelRMSD)
					{
						if (majorAtom == "N")
						{
//							continue;
						}

						addDampening(bond, 0.1, 0.1);
						addTorsionBlur(bond, 0.1, 0.1);
					}

					setJobName("magic_axis_" + i_to_str(resNum) + "_" + bond->shortDesc());
					setSilent();
					addSampledBackbone(getPolymer(), resNum, resNum + 20);
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
				addSampledBackbone(getPolymer(), resNum, resNum + 3);

				setScoreType(ScoreTypeModelPos);
				sample();

				std::cout << "Res " << resNum << " bond " << bond->shortDesc();
				std::cout << " : bFactor " << bond->getMeanSquareDeviation();
				std::cout << std::endl;
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
					setCycles(10);
					setupTorsionSet(bond, k, 5, resNum, 0.05, 0.02, true);
					setCrystal(target);
					sample();
				}
/*
				setupNelderMead();
				setupTorsionSet(bond, k, 4, resNum, 1.0, 0.05);
				setCrystal(target);
				sample();
*/
				setupNelderMead();
				setupTorsionSet(bond, k, 4, resNum, 0.3, 0.02);
//				setupTorsionSet(bond, k, 3, resNum, 0.05, 0.02, true);
				setCrystal(target);
				sample();
			}


			bond->setActiveGroup(0);

		}

	}
}

AtomPtr Backbone::betaCarbonTorsionAtom()
{
	/* What is the major atom of the bond describing CA? */

	AtomPtr ca = findAtom("CA");
	ModelPtr model = ca->getModel();

	if (model->getClassName() == "Bond")
	{
		BondPtr bond = ToBondPtr(model);
		return bond->getMajor();
	}

	return AtomPtr();
}
