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

	const char *atoms[] = {"CB", "CG", "CG1", "CG2", "CD", "CE", "NZ", "OG1"};
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

			if (!bond->isNotJustForHydrogens())
			{
				continue;
			}

			std::string preAtom = bond->getMajor()->getAtomName();
			int groups = bond->downstreamAtomGroupCount();

			if (preAtom != "CA" && bond->isUsingTorsion()
				&& (rType == RefinementFine))
			{

				for (int k = 0; k < groups; k++)
				{
					bond->setActiveGroup(k);

					setupNelderMead();
					addTorsionNextBlur(bond, 0.1, deg2rad(0.1));

					for (int j = 0; j < bond->downstreamAtomCount(k); j++)
					{
						addSampled(bond->downstreamAtom(k, j));
					}

					setJobName("compensate_" + preAtom + "_"
							   + atom + "_" + i_to_str(resNum));
					setCrystal(target);
					sample();
				}

				bond->setActiveGroup(0);
			}

			ModelPtr preModel = bond->getParentModel();

			if (preModel->getClassName() == "Bond" && (rType == RefinementFine))
			{
				setupNelderMead();
				reportInDegrees();

				BondPtr preBond = std::static_pointer_cast<Bond>(preModel);

				addBendAngle(bond, deg2rad(0.2), deg2rad(0.1));
				addBendBlur(bond, deg2rad(0.05), deg2rad(0.01));

				addSampled(bond->getMinor());

				setJobName("bend_" + preAtom + "_" + atom + "_" + i_to_str(resNum));
				setCrystal(target);
				sample();
			}

			if (bond->isUsingTorsion())
			{
				for (int k = 0; k < groups; k++)
				{
					bond->setActiveGroup(k);

					if (rType == RefinementFine)
					{
						setupNelderMead();
					}
					else
					{
						setupGrid();
					}

					reportInDegrees();

					if (rType == RefinementFine)
					{
						addTorsion(bond, deg2rad(1.0), deg2rad(0.2));
						addTorsionBlur(bond, deg2rad(1.0), deg2rad(0.2));
					}
					else
					{
						addTorsion(bond, deg2rad(150), deg2rad(1.0));
					}

					for (int j = 0; j < bond->downstreamAtomCount(k); j++)
					{
						addSampled(bond->downstreamAtom(k, j));
					}

					setJobName("torsion_" + preAtom + "_" + atom + "_g" +
							   i_to_str(k) + "_" + i_to_str(resNum));
					setCrystal(target);
					sample();
				}
/*
				if (groups > 1 && rType == RefinementFine)
				{
					for (int k = 0; k < groups - 1; k++)
					{
						bond->setActiveGroup(k);
						setupGrid();
						setJointSampling();

						Bond::setOccupancy(&*bond, 0.5);
						addOccupancy(bond, 1.0, 0.01);

						for (int j = 0; j < bond->downstreamAtomGroupCount(); j++)
						{
							for (int l = 0; l < bond->downstreamAtomCount(j); l++)
							{
								addSampled(bond->downstreamAtom(j, l));
							}
						}

						setJobName("occupancy_" + preAtom + "_" + atom + "_" + i_to_str(resNum));
						setCrystal(target);
						sample();
					}
				}
*/
				bond->setActiveGroup(0);

			}

			if (preAtom == "CB" && myAtoms.size() <= 1)
			{
				bond->splitBond();
//				int newGroup = bond->downstreamAtomGroupCount() - 1;

				refine(target, RefinementBroad);
			}
		}

	}
}


/*
 if (preAtom == "CB")
 {
 ModelPtr preModel = bond->getParentModel();
 if (preModel->getClassName() != "Bond")
 {
 break;
 }

 BondPtr preBond = std::static_pointer_cast<Bond>(preModel);

 setupGrid();
 addTorsion(preBond, deg2rad(180.0), deg2rad(4.0));
 addBendAngle(bond, deg2rad(90.0), deg2rad(1.0));

 addSampled(bond->getMinor());

 setJobName("grid_" + preAtom + "_" + atom + "_" + i_to_str(resNum));
 setCrystal(target);
 sample();
 }
 */