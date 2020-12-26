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
	flexer->_workspace.tBonds->setFine(true);
	flexer->_workspace.tMap->setFine(true);
	flexer->_workspace.tScore->setFine(true);

	flexer->_workspace.tBonds->start();
	flexer->_atomGroup->refreshPositions();
	flexer->_workspace.tBonds->stop();
	
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
	_mapList.clear();
	_paramList.clear();
}

void FlexGlobal::prepareComparisons(RefinementStrategyPtr str)
{
	/*
	prepareWorkspace();
	MapScoreWorkspace *ws = &getWorkspace();
	FlexGlobal::score(this);
	VagFFTPtr s = VagFFTPtr(new VagFFT(*ws->segment));
	std::map<int, VagFFTPtr> maps;
	_paramCorrel.clear();
	_paramList.clear();

	for (size_t i = 0; i < str->parameterCount(); i++)
	{
		Parameter p = str->getParamObject(i);
		VagFFTPtr diff = VagFFTPtr(new VagFFT(*ws->segment));
		diff->multiplyFinal(-1);
		double orig = p.start_value;
		double newval = orig + p.step_size / 10.;
		(p.setter)(p.object, newval); 
		FlexGlobal::score(this);
		diff->addSimple(ws->segment);
		(p.setter)(p.object, orig); 
		FlexGlobal::score(this);
		maps[i] = diff;
		_paramList.push_back(p);
	}

	for (size_t i = 0; i < str->parameterCount(); i++)
	{
		Parameter p1 = str->getParamObject(i);

		for (size_t j = 0; j < str->parameterCount(); j++)
		{
			if (i == j)
			{
				continue;
			}

			Parameter p2 = str->getParamObject(j);

			VagFFTPtr diff1 = VagFFTPtr(new VagFFT(*(ws->segment)));
			diff1->multiplyFinal(-1);
			VagFFTPtr diff2 = VagFFTPtr(new VagFFT(*diff1));
			void *obj1 = p1.object;
			void *obj2 = p2.object;
			double start1 = p1.start_value;
			double start2 = p2.start_value;
			double inc1 = p1.step_size / 10.;
			double inc2 = p2.step_size / 10.;

			(*p1.setter)(obj1, start1 + inc1); 

			FlexGlobal::score(this);
			diff1->addSimple(ws->segment);

			(*p1.setter)(obj1, start1); 
			(*p2.setter)(obj2, start2 + inc2); 

			FlexGlobal::score(this);
			diff2->addSimple(ws->segment);

			CorrelData cd  = empty_CD();
			CorrelData cd2 = empty_CD();
			for (size_t k = 0; k < maps[i]->nn(); k++)
			{
				add_to_CD(&cd2, diff1->getReal(k), diff2->getReal(k));
				add_to_CD(&cd, maps[i]->getReal(k), maps[j]->getReal(k));
			}

			(*p2.setter)(obj2, start2); 
			FlexGlobal::score(this);

			double correl = evaluate_CD(cd);
			double correl2 = evaluate_CD(cd2);
			std::cout << correl << " " << correl2 <<  std::endl;
			_paramCorrel[i][j] = correl2;

			CorrelData cd = empty_CD();

			for (size_t k = 0; k < maps[i]->nn(); k++)
			{
				double r1 = maps[i]->getReal(k);
				double r2 = maps[j]->getReal(k);
				add_to_CD(&cd, r1, r2);
			}

			double correl = evaluate_CD(cd);
			std::cout << correl << std::endl;
			_paramCorrel[i][j] = correl;
		}
	}
			*/
}

bool equal_params(Parameter &p1, Parameter &p2)
{
	return (p1.object == p2.object && p1.getter == p2.getter);
}

double FlexGlobal::compareParams(void *object, Parameter &p1, Parameter &p2)
{
	FlexGlobal *me = static_cast<FlexGlobal *>(object);

	int i1 = -1; int i2 = -1;
	for (size_t i = 0; i < me->_paramList.size(); i++)
	{
		if (equal_params(me->_paramList[i], p1))
		{
			i1 = i;
		}
		else if (equal_params(me->_paramList[i], p2))
		{
			i2 = i;
		}
	}

	MapScoreWorkspace *ws = &me->getWorkspace();
	VagFFTPtr diff1;
	VagFFTPtr diff2;

	void *obj1 = p1.object;
	void *obj2 = p2.object;
	double start1 = p1.start_value;
	double start2 = p2.start_value;
	double inc1 = p1.step_size / 10.;
	double inc2 = p2.step_size / 10.;

	if (i1 < 0)
	{
		diff1 = VagFFTPtr(new VagFFT(*(ws->segment)));
		diff1->multiplyFinal(-1);
		(*p1.setter)(obj1, start1 + inc1); 
		FlexGlobal::score(me);
		diff1->addSimple(ws->segment);
		(*p1.setter)(obj1, start1); 
		FlexGlobal::score(me);

		me->_paramList.push_back(p1);
		me->_mapList.push_back(diff1);
		i1 = me->_mapList.size() - 1;
	}

	if (i2 < 0)
	{
		diff2 = VagFFTPtr(new VagFFT(*ws->segment));
		diff2->multiplyFinal(-1);
		(*p2.setter)(obj2, start2 + inc2); 
		FlexGlobal::score(me);
		diff2->addSimple(ws->segment);
		(*p2.setter)(obj2, start2); 
		FlexGlobal::score(me);

		me->_paramList.push_back(p2);
		me->_mapList.push_back(diff2);
		i2 = me->_mapList.size() - 1;
	}

	diff1 = me->_mapList[i1];
	diff2 = me->_mapList[i2];

	CorrelData cd = empty_CD();
	for (size_t i = 0; i < diff1->nn(); i++)
	{
		add_to_CD(&cd, diff1->getReal(i), diff2->getReal(i));
	}
	
	double correl = evaluate_CD(cd);
//	std::cout << correl << std::endl;
	return correl;
}

void FlexGlobal::reportTimings()
{
	std::cout << std::endl;
	_workspace.tBonds->report();
	_workspace.tMap->report();
	_workspace.tScore->report();

}
