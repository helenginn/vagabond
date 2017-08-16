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
	int resNum = getMonomer()->getResidueNum();

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
			BondPtr bond = std::static_pointer_cast<Bond>(myAtom->getModel());
			std::string majorAtom = bond->getMajor()->getAtomName();

			if (_resNum == 120 && majorAtom == "CG" && myAtoms.size() <= 1)
			{
				bond->splitBond();
				ModelPtr upModel = bond->getParentModel();
				BondPtr upBond = std::static_pointer_cast<Bond>(upModel);
				upBond->setActiveGroup(1);

				setupGrid();
				reportInDegrees();

				addTorsion(upBond, deg2rad(270), deg2rad(1.0));

				for (int j = 0; j < upBond->downstreamAtomCount(1); j++)
				{
					addSampled(upBond->downstreamAtom(1, j));
				}

				for (int j = 0; j < upBond->extraTorsionSampleCount(1); j++)
				{
					addSampled(upBond->extraTorsionSample(1, j));
				}

				setJobName("torsion_" + majorAtom + "_" + atom + "_g" +
						   i_to_str(1) + "_" + i_to_str(resNum));
				setCrystal(target);
				sample();

				upBond->setActiveGroup(0);
			}

			if (!bond->isNotJustForHydrogens() || bond->isFixed())
			{
				continue;
			}

			int groups = bond->downstreamAtomGroupCount();

			if (bond->isUsingTorsion() && (rType == RefinementFine))
			{
				for (int k = 0; k < groups; k++)
				{
					bond->setActiveGroup(k);
					setupGrid();
					Bond::setTorsionNextBlur(&*bond, -0.5);
					addTorsionNextBlur(bond, 1.0, 0.02);

					for (int j = 0; j < bond->downstreamAtomCount(k); j++)
					{
						addSampled(bond->downstreamAtom(k, j));
					}

					setJobName("compensate_" + majorAtom + "_"
							   + atom + "_" + i_to_str(resNum));
					setCrystal(target);
					sample();
				}

				bond->setActiveGroup(0);
			}

			ModelPtr preModel = bond->getParentModel();

			if (false && preModel->getClassName() == "Bond" && (rType == RefinementFine))
			{
				setupNelderMead();
				reportInDegrees();

				BondPtr preBond = std::static_pointer_cast<Bond>(preModel);

		//		addBendAngle(bond, deg2rad(0.2), deg2rad(0.1));
		//		addBendBlur(bond, deg2rad(0.05), deg2rad(0.01));

				addSampled(bond->getMinor());

				setJobName("bend_" + majorAtom + "_" + atom + "_" + i_to_str(resNum));
				setCrystal(target);
				sample();
			}

			if (bond->isUsingTorsion() && bond->isNotJustForHydrogens())
			{
				if (rType == RefinementBroad)
				{
					for (int k = 0; k < groups; k++)
					{
						setupDoubleTorsion(bond, k, 1, resNum, 360, 8);
						setCrystal(target);
						sample();

						setupDoubleTorsion(bond, k, 1, resNum, 30, 2);
						setCrystal(target);
						sample();
					}
				}
				else
				{
					for (int k = 0; k < groups; k++)
					{
						bond->setActiveGroup(k);

						setupGrid();
						reportInDegrees();

						addTorsion(bond, deg2rad(0.5), deg2rad(0.2));
						addTorsionBlur(bond, deg2rad(0.2), deg2rad(0.1));

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
			}
		}

	}
}
