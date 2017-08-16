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
#include "FileReader.h"

void Backbone::refine(CrystalPtr target, RefinementType rType)
{
	if (rType == RefinementBroad)
	{
	//	return;
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

			if (model->getClassName() != "Bond")
			{
				continue;
			}

			BondPtr bond = std::static_pointer_cast<Bond>(model);

			if (!bond->isNotJustForHydrogens() || bond->isFixed())
			{
				continue;
			}

			std::string majorAtom = bond->getMajor()->getAtomName();
			int groups = bond->downstreamAtomGroupCount();

			if (bond->isUsingTorsion() && (rType == RefinementFine))
			{
				for (int k = 0; k < groups; k++)
				{
					bond->setActiveGroup(k);
					setupGrid();
					Bond::setTorsionNextBlur(&*bond, 0.5);
					addTorsionNextBlur(bond, 1.0, 0.02);
					addSampled(bond->getMinor());

					for (int j = 0; j < bond->downstreamAtomCount(k); j++)
					{
						addSampled(bond->downstreamAtom(k, j));
					}

					setJobName("compensate_" + majorAtom + "_"
							   + atom + "_" + i_to_str(resNum));
					setCrystal(target);
					sample();
				}

//				bond->setActiveGroup(0);
			}

			ModelPtr preModel = bond->getParentModel();
/*
			if (preModel->getClassName() == "Bond" && (rType == RefinementFine))
			{
				setupNelderMead();
				reportInDegrees();

				BondPtr preBond = std::static_pointer_cast<Bond>(preModel);

				addBendAngle(bond, deg2rad(0.2), deg2rad(0.1));
				addBendBlur(bond, deg2rad(0.05), deg2rad(0.01));

				addSampled(bond->getMinor());

				setJobName("bend_" + majorAtom + "_" + atom + "_" + i_to_str(resNum));
				setCrystal(target);
				sample();
			}
*/
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
							addTorsion(bond, deg2rad(0.5), deg2rad(0.2));
							addTorsionBlur(bond, deg2rad(0.2), deg2rad(0.2));
						}
						else
						{
							addTorsion(bond, deg2rad(360), deg2rad(1.0));
						}

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

			continue;

			if (_resNum == 123 && majorAtom == "N" && myAtoms.size() <= 1)
			{
				bond->splitBond();

				refine(target, RefinementBroad);
			}
		}
		
	}
}
