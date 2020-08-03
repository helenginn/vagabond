//
//    RefinementLBFGS.h
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

#ifndef __vagabond__RefinementLBFGS__
#define __vagabond__RefinementLBFGS__

#include "lbfgs.h"
#include "RefinementStrategy.h"

typedef std::vector<lbfgsfloatval_t> LbfgsVector;

class RefinementLBFGS : public RefinementStrategy
{
public:
	RefinementLBFGS();

	void setGradientRefresh(void *gradObj, Getter func)
	{
		_gradObj = gradObj;
		_func = func;
	}
	
	virtual void refine();
private:
	bool hasAllGradients();

	static double evaluate(void *instance,
	                       const lbfgsfloatval_t *x,
	                       lbfgsfloatval_t *g,
	                       const int n,
	                       const lbfgsfloatval_t step);

	static int progress( void *instance, const lbfgsfloatval_t *x,
	                    const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
	                    const lbfgsfloatval_t xnorm, 
	                    const lbfgsfloatval_t gnorm,
	                    const lbfgsfloatval_t step, int n, int k, int ls);

	void copyInGradientValues(lbfgsfloatval_t *g);
	void copyInStartValues();
	void copyOutValues(const lbfgsfloatval_t *x);

	void *_gradObj;
	Getter _func;
	lbfgsfloatval_t _fx;
	LbfgsVector _xs;
	LbfgsVector _gs;
};

#endif
