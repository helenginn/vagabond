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

double FlexGlobal::localScaleParam(Parameter &p1)
{
	double start1 = (*p1.getter)(p1.object);
	double change1 = p1.step_size + start1;

	score(this);
	FFTPtr groundState = FFTPtr(new FFT(*(_workspace.segment)));

	(*p1.setter)(p1.object, change1);
	score(this);
	(*p1.setter)(p1.object, start1);

	FFTPtr p1State = FFTPtr(new FFT(*(_workspace.segment)));
	FFT::addSimple(p1State, groundState, -1);
	
	double sum_all = groundState->sumReal();
	double sum_frac = p1State->sumReal();
	
	if (sum_frac <= 1e-6)
	{
		p1.step_size = 0;
		return 0;
	}
	
	double frac = sum_frac / sum_all;
	
	double scale = 0.05 / frac;
	p1.step_size *= scale;
	
	return scale;
}

double FlexGlobal::localCompareParams(Parameter &p1, Parameter &p2)
{
	double start1 = (*p1.getter)(p1.object);
	double start2 = (*p2.getter)(p2.object);

	double change1 = p1.step_size + start1;
	double change2 = p2.step_size + start2;

	score(this);
	
	FFTPtr groundState = FFTPtr(new FFT(*(_workspace.segment)));

	(*p1.setter)(p1.object, change1);
	score(this);
	FFTPtr p1State = FFTPtr(new FFT(*(_workspace.segment)));
	(*p1.setter)(p1.object, start1);
	FFT::addSimple(p1State, groundState, -1);

	(*p2.setter)(p2.object, change2);
	score(this);
	(*p2.setter)(p2.object, start2);

	FFTPtr p2State = FFTPtr(new FFT((*_workspace.segment)));
	FFT::addSimple(p2State, groundState, -1);

	std::vector<double> xs, ys;
	for (int i = 0; i < p1State->nn; i++)
	{
		xs.push_back(p1State->data[i][0]);
		ys.push_back(p2State->data[i][0]);
	}
	
	double correl = correlation(xs, ys);
	
	if (correl != correl)
	{
		correl = 0;
	}

	return correl;
}

double FlexGlobal::compareParams(void *obj, Parameter &p1, Parameter &p2)
{
	FlexGlobal *me = static_cast<FlexGlobal *>(obj);
	double score = me->localCompareParams(p1, p2);
	
	return score;
}

double FlexGlobal::scaleParam(void *obj, Parameter &p1)
{
	FlexGlobal *me = static_cast<FlexGlobal *>(obj);
	double score = me->localScaleParam(p1);
	
	return score;
}
