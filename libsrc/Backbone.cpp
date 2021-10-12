//
//  Backbone.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Backbone.h"
#include "Monomer.h"
#include "Options.h"
#include "Bond.h"
#include "Atom.h"
#include "Absolute.h"
#include <hcsrc/FileReader.h>
#include "Sidechain.h"

bool Backbone::shouldRefineAtom(AtomPtr atom)
{
	return true;
}

void Backbone::refine(CrystalPtr target, RefinementType rType)
{
	if (!isFullyTied())
	{
		return;
	}

	if (rType == RefinementSidechain || rType == RefinementSidePos)
	{
		return;
	}

	*_stream << getMonomer()->getResCode() << std::flush;
	
	if (_timesRefined == 0)
	{
		_excludeO = true;
	}
	
	if (!paramCount())
	{
		double range = 2.;

		if (_timesRefined >= 6)
		{
			range = 0.5;
		}

		switch (rType)
		{
			case RefinementModelPos:
			addParamType(ParamOptionTorsion, range);
			addParamType(ParamOptionMaxTries, 60);
			addParamType(ParamOptionNumBonds, 3);
			break;
			
			case RefinementSavedPos:
			addParamType(ParamOptionTorsion, range);
			addParamType(ParamOptionBondAngle, range / 1);
			addParamType(ParamOptionNumBonds, 3);
			break;
			
			case RefinementCrude:
			addParamType(ParamOptionTorsion, 2.0);
			addParamType(ParamOptionTopLevelOnly, 2.0);
			addParamType(ParamOptionFirstOccupancy, 2.0);
			addParamType(ParamOptionSVD, 1);
			addParamType(ParamOptionCycles, 120);
			addParamType(ParamOptionNumBonds, 10);
			addParamType(ParamOptionExtraAtoms, 2);
			break;
			
			case RefinementFine:
			addParamType(ParamOptionTorsion, 12.0);
			addParamType(ParamOptionTopLevelOnly, 2.0);
			addParamType(ParamOptionFirstOccupancy, 2.0);
			addParamType(ParamOptionCycles, 20);
			addParamType(ParamOptionNumBonds, 2);
			addParamType(ParamOptionMaxTries, 5);
			break;

			default:
			break;
		}
		
		if (_timesRefined >= 5)
		{
			switch (rType)
			{
				case RefinementModelPos:
				addParamType(ParamOptionNumBonds, 8);
				addParamType(ParamOptionCycles, 100);
				addParamType(ParamOptionMaxTries, 15);
				addParamType(ParamOptionThorough, 1);
				addParamType(ParamOptionSVD, 1);
				break;

				case RefinementFine:
//				addParamType(ParamOptionBondAngle, 0.1);
				addParamType(ParamOptionCirclePortion, 0.1);

				default:
				break;
			}
		}
		
		if (_timesRefined == 8)
		{
			AtomPtr c = findAtom("C");
			AtomPtr o = findAtom("O");
			if (c && o)
			{
				double disp = o->posDisplacement();
				BondPtr b = ToBondPtr(o->getModel());
				b->flipPyramid();
				o->getExplicitModel()->getFinalPositions();
				double new_disp = o->posDisplacement();
				bool change = new_disp < disp;

				if (new_disp > disp)
				{
					b->flipPyramid();
				}

				_excludeO = false;
			}
		}
	}

	MonomerPtr monomer = getMonomer();
	_includeForRefine = monomer->includingInRefinement();

	AtomGroup::refine(target, rType);
	clearParams();
	clearIncludeForRefinements();
}

AtomPtr Backbone::betaCarbonTorsionAtom()
{
	/* What is the major atom of the bond describing CA? */

	AtomPtr ca = findAtom("CA");

	if (!ca)
	{
		return AtomPtr();
	}

	ModelPtr model = ca->getModel();

	if (model->getClassName() == "Bond")
	{
		BondPtr bond = ToBondPtr(model);
		return bond->getMajor();
	}

	return AtomPtr();
}

void Backbone::addProperties()
{
	addIntProperty("res_num", &_resNum);
	AtomGroup::addProperties();
}
