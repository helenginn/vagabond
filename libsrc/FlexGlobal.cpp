//
//  FlexGlobal.cpp
//  vagabond
//
//  Created by Helen Ginn on 27/12/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#include "FlexGlobal.h"
#include "AtomGroup.h"
#include "Atom.h"
#include "Bond.h"
#include "maths.h"
#include "Anisotropicator.h"
#include "Shouter.h"
#include "FFT.h"
#include "Crystal.h"

FlexGlobal::FlexGlobal()
{
	_prepared = false;
	_recip = false;
}

void FlexGlobal::prepareWorkspace()
{
	if (_prepared) return;

	if (!_crystal)
	{
		shout_at_helen("Helen trying to fit electron density"\
		               "to missing crystal");	
	}

	setup_space(&_workspace);

	_workspace.crystal = _crystal;
	_workspace.selectAtoms = _atomGroup;

	_workspace.selectAtoms->addParamType(ParamOptionStep, 2);
	
	_prepared = true;
}

double FlexGlobal::score(void *object)
{
	FlexGlobal *flexer = static_cast<FlexGlobal *>(object);
	flexer->prepareWorkspace();

	flexer->_atomGroup->refreshPositions();
	
	double score = 0;
	if (flexer->_recip)
	{
		score = AtomGroup::scoreWithReciprocal(&flexer->_workspace);
	}
	else
	{
		score = AtomGroup::scoreWithMapGeneral(&flexer->_workspace);
	}

	return score;
}

void FlexGlobal::plot(std::string filename)
{
	_workspace.filename = filename;
	AtomGroup::scoreWithMapGeneral(&_workspace, true);
}

void FlexGlobal::recalculateConstant()
{
	_workspace.recalc = true;
}

