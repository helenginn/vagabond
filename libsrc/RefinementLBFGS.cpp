//
//    RefinementLBFGS.cpp
//    Vagabond - stuff

//    Copyright (C) 2018  Helen Ginn
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "RefinementLBFGS.h"
#include <iostream>
#include <iomanip>

RefinementLBFGS::RefinementLBFGS()
{
	_fx = 0;
	_func = NULL;
}

int RefinementLBFGS::progress( void *instance, const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step, int n, int k, int ls)
{
    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}

double RefinementLBFGS::evaluate(void *instance,
                                 const lbfgsfloatval_t *x,
                                 lbfgsfloatval_t *g,
                                 const int n,
                                 const lbfgsfloatval_t step)
{
	RefinementLBFGS *me = static_cast<RefinementLBFGS *>(instance);

	// put the values of x into the setters.
	me->copyOutValues(x);
	
	if (me->_func)
	{
		(*me->_func)(me->_gradObj);
	}
	
	// compute new values of gradients
	me->copyInGradientValues(g);

	// return a new fx evaluation.
	double eval = me->evaluationFunction(me->evaluateObject);
	std::cout << "Evaluation: " << std::setprecision(10) << eval << std::endl;
	
	/*
	std::cout << "Evaluation list:" << std::endl;
	
	for (int i = 0; i < n; i++)
	{
		std::cout << std::setprecision(10) <<
		x[i] << ", " << g[i] << std::endl;
	}
	*/
	
	
	return eval;
}

void RefinementLBFGS::copyOutValues(const lbfgsfloatval_t *x)
{
	for (int i = 0; i < parameterCount(); i++)
	{
		setValueForParam(i, x[i]);
	}
}

void RefinementLBFGS::copyInStartValues()
{
	for (int i = 0; i < parameterCount(); i++)
	{
		_xs[i] = getValueForParam(i);
	}
}

void RefinementLBFGS::copyInGradientValues(lbfgsfloatval_t *g)
{
	for (int i = 0; i < parameterCount(); i++)
	{
		g[i] = getGradientForParam(i);
	}
}

bool RefinementLBFGS::hasAllGradients()
{
	for (int i = 0; i < parameterCount(); i++)
	{
		if (_params[i].gradient == NULL)
		{
			return false;
		}
	}

	return true;
}

void RefinementLBFGS::refine()
{
	RefinementStrategy::refine();
	lbfgs_parameter_t param;
	
	if (!hasAllGradients())
	{
//		std::cout << "LBFGS needs gradients atm." << std::endl;
//		return;
	}
	
	std::cout << "Starting LBFGS" << std::endl;

	lbfgs_parameter_init(&param);
	param.epsilon = 1e-10;
	param.max_iterations = 10;
	
	int count = parameterCount();
	
	_xs.resize(count);
	_gs.resize(count);
	
	copyInStartValues();
//	copyInGradientValues(&_gs[0]);

	if (_func == NULL)
	{
		std::cout << "Yelp" << std::endl;
	}
	
	int result = lbfgs(count, &_xs[0], &_fx, 
	                   evaluate, progress, this, &param);

	finish();
}
