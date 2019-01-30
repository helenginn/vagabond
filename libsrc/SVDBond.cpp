// Vagabond
// Copyright (C) 2017-2018 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#include <iostream>
#include <string.h>
#include "SVDBond.h"
#include "../libica/svdcmp.h"
#include "Atom.h"
#include "RefinementStrategy.h"
#include "Bond.h"
#include "Whack.h"

SVDBond::SVDBond(BondEffects &effects, std::vector<BondPtr> &bonds,
                 std::vector<AtomPtr> &atoms)
{
	_effects = effects;
	_bonds = bonds;
	_atoms = atoms;
	_svd = NULL;
}

void SVDBond::performSVD(BondBondCC *ccs)
{
	prepareMatrix(&_svd);
	prepareMatrix(&_original);
	prepareMatrix(&_v);
	prepareVector(&_w);
	
	if (ccs)
	{
		for (int i = 0; i < _bonds.size(); i++)
		{
			BondPtr bi = _bonds[i];

			for (int j = 0; j < _bonds.size(); j++)
			{
				BondPtr bj = _bonds[j]; // giggle
				double correl = (*ccs)[bi][bj];
				
				if (i == j)
				{
					correl = 0;
				}
				
				_svd[i][j] = correl;
			}
		}
	}
	else
	{
		populateMatrix();
	}
	
	copyMatrix(_svd, _original);
	svdMagic();
	report();
}

void SVDBond::prepareVector(double **ptr)
{
	size_t svd_dims = _bonds.size();

	*ptr = (double *)malloc(sizeof(double) * svd_dims);
}

void SVDBond::prepareMatrix(double ***ptr)
{
	size_t svd_dims = _bonds.size();

	*ptr = (double **)malloc(sizeof(double *) * svd_dims);
	for (int i = 0; i < svd_dims; i++)
	{
		(*ptr)[i] = (double *)malloc(sizeof(double) * svd_dims);
	}
}

void SVDBond::copyMatrix(double **from, double **to)
{
	size_t svd_dims = _bonds.size();

	for (int i = 0; i < svd_dims; i++)
	{
		memcpy(to[i], from[i], sizeof(double) * svd_dims);
	}
}

void SVDBond::populateMatrix()
{
	for (int i = 0; i < _bonds.size(); i++)
	{
		BondPtr bi = _bonds[i];
		AtomTarget ai = _effects[bi];

		for (int j = 0; j < _bonds.size(); j++)
		{
			BondPtr bj = _bonds[j];
			AtomTarget aj = _effects[bj];
			
			double total = 0;
			for (int k = 0; k < _atoms.size(); k++)
			{
				double add = ai[_atoms[k]] * aj[_atoms[k]];
				total += add;
			}
			
			_svd[i][j] = total;
		}
	}
}

void SVDBond::svdMagic()
{
	size_t dims = _bonds.size();
	int success = svdcmp((mat)_svd, dims, dims, (vect)_w, (mat)_v);
	
	if (!success)
	{
		return;
	}
}

void SVDBond::report()
{
	double cuml = 0;
	for (int i = 0; i < _bonds.size(); i++)
	{
		cuml += _w[i];
	}
	
	_wTotal = cuml; 
	
	double pc99 = cuml * 0.50;
	cuml = 0;
	_num = 0;
	for (int i = 0; i < _bonds.size(); i++)
	{
		cuml += _w[i];
		
		if (cuml > pc99)
		{
			_num = i;
			break;
		}
	}
	
	if (_num < 5)
	{
		num = 5;
		
		if (_bonds.size() < 5)
		{
			num = _bonds.size();
		}
	}
	
	/* Set up parameters for the top clusters from SVD */
	for (int i = 0; i < _num; i++)
	{
		ParamSVDSet set;
		Param *pw = new Param();
		Param::setValue(pw, 0);

		Param *pk = new Param();
		Param::setValue(pk, 0);
		
		set.pWhack = pw;
		set.pKick = pk;
		set.rowPtr = _svd[i];

		_params.push_back(set);
	}
	
	/* Get the baseline values for kicks and whacks */

	for (int i = 0; i < _bonds.size(); i++)
	{
		BondPtr bond = _bonds[i];
		WhackPtr whack = bond->getWhack();
		if (!whack) continue;
		
		double w = Whack::getWhack(&*whack);
		double k = Whack::getKick(&*whack);
		
		BondParamPair pair;
		pair.whack = w;
		pair.kick = k;
		_bondBase[bond] = pair;
	}
}

void SVDBond::addToStrategy(RefinementStrategyPtr strategy)
{
	double inv = 0.5 / _wTotal;
	double tol = inv / 100;
	
	for (int i = 0; i < _params.size(); i++)
	{
		strategy->addParameter(_params[i].pWhack, Param::getValue,
		                       Param::setValue, inv, tol);
		strategy->addParameter(_params[i].pKick, Param::getValue,
		                       Param::setValue, inv, tol);
	}
}

void SVDBond::applyParameters()
{
	/* for every bond... */
	for (int i = 0; i < _bonds.size(); i++)
	{
		BondPtr bi = _bonds[i];
		
		if (_bondBase.count(bi) == 0)
		{
			continue;
		}
		
		double w = _bondBase[bi].whack;
		double k = _bondBase[bi].kick;
		
//		std::cout << "base: " << w << ", " << k << std::endl;
		
		/* to store the total contributions to add to whack/kick */
		double w_more = 0;
		double k_more = 0;

		/* loop round each cluster... */
		for (int j = 0; j < _params.size(); j++)
		{
			/* get the relative portion of this bond's identification
			 * with this cluster */
			
			double total = 0;
			for (int k = 0; k < _bonds.size(); k++)
			{
				BondPtr bk = _bonds[k];
				double add = _original[i][k] * _params[j].rowPtr[k];
				total += add;
			}
			
			/* get value of the whack/kick */
			double wVal = Param::getValue(_params[j].pWhack);
			double kVal = Param::getValue(_params[j].pKick);

			wVal *= total;
			kVal *= total;
			
			w_more += wVal;
			k_more += kVal;
		}
		
		w += w_more;
		k += k_more;

//		std::cout << "now: " << w << ", " << k << std::endl;
		
		WhackPtr whack = bi->getWhack();
		Whack::setWhack(&*whack, w);
		Whack::setKick(&*whack, k);
	}

}

void SVDBond::cleanupSVD(double ***ptr)
{
	if (*ptr == NULL)
	{
		return;
	}
	
	size_t svd_dims = _bonds.size();

	for (int i = 0; i < svd_dims; i++)
	{
		free((*ptr)[i]);
		(*ptr)[i] = NULL;
	}
	
	free(*ptr);
	*ptr = NULL;
}

SVDBond::~SVDBond()
{
	cleanupSVD(&_svd);
	cleanupSVD(&_v);
	
	free(_w); _w = NULL;
	
	for (int i = 0; i < _params.size(); i++)
	{
		delete _params[i].pWhack;
		delete _params[i].pKick;
	}
}
