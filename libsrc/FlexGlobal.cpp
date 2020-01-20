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
#include "Anisotropicator.h"
#include "Shouter.h"

FlexGlobal::FlexGlobal()
{
	_prepared = false;
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
	
	_prepared = true;
	AtomGroup::scoreWithMapGeneral(&_workspace, false);
}

double FlexGlobal::score(void *object)
{
	FlexGlobal *flexer = static_cast<FlexGlobal *>(object);
	flexer->prepareWorkspace();

	flexer->_atomGroup->propagateChange();
	double score = AtomGroup::scoreWithMapGeneral(&flexer->_workspace);

	return score;
}

void FlexGlobal::plot(std::string filename)
{
	_workspace.filename = filename;
	AtomGroup::scoreWithMapGeneral(&_workspace, true);
}

double FlexGlobal::localCompareParams(Parameter &p1, Parameter &p2)
{
	return 0;
}

double FlexGlobal::compareParams(void *obj, Parameter &p1, Parameter &p2)
{
	FlexGlobal *me = static_cast<FlexGlobal *>(obj);
	double score = me->localCompareParams(p1, p2);
	
	return score;
}
