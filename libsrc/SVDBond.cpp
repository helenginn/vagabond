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
#include "Polymer.h"
#include "Anchor.h"
#include "Whack.h"
#include "CSV.h"

SVDBond::SVDBond(std::vector<BondPtr> &bonds, std::vector<AtomPtr> &atoms)
{
	_silent = true;
	_doTorsion = false;
	_bonds = bonds;
	_atoms = atoms;
	_svd = NULL;
}

SVDBond::SVDBond(std::vector<AtomPtr> &atoms)
{
	_silent = true;
	_doTorsion = false;
	_atoms = atoms;
	_svd = NULL;
}

vec3 angle_effect_on_pos(vec3 atom_pos, mat3x3 &bond_basis, vec3 &bond_pos)
{
	mat3x3 trans = mat3x3_transpose(bond_basis);
	/* Remove translation component due to bond's location */
	vec3_subtract_from_vec3(&atom_pos, bond_pos);
	
	/* Multiply position by transpose of bond basis to get location in
	 * 'identity' coordinates */
	mat3x3_mult_vec(trans, &atom_pos);
	
	/* Set the component in the direction of the bond to zero */
	atom_pos.z = 0;
	
	/* Multiply this by bond basis to get the important vector into
	 * crystal coordinates again */
	mat3x3_mult_vec(bond_basis, &atom_pos);
	
	return atom_pos;
}

vec3 bond_effect_on_pos(vec3 atom_pos, mat3x3 &bond_basis, vec3 &bond_pos)
{
	mat3x3 trans = mat3x3_transpose(bond_basis);
	/* Remove translation component due to bond's location */
	vec3_subtract_from_vec3(&atom_pos, bond_pos);
	
	/* Multiply position by transpose of bond basis to get location in
	 * 'identity' coordinates */
	mat3x3_mult_vec(trans, &atom_pos);
	
	/* Set the component in the direction of the bond to zero */
	atom_pos.z = 0;
	
	/* Multiply this by bond basis to get the important vector into
	 * crystal coordinates again */
	mat3x3_mult_vec(bond_basis, &atom_pos);
	
	return atom_pos;
}

double SVDBond::compareKicks(BondPtr a, BondPtr b)
{
	BondPtr aChild = a->downstreamBond(0, 0);
	BondPtr bChild = b->downstreamBond(0, 0);
	std::vector<BondSample> as = aChild->getFinalPositions();
	std::vector<BondSample> bs = bChild->getFinalPositions();

	double total = 0;

	for (size_t i = 0; i < as.size(); i++)
	{
		double mult = as[i].kickValue * bs[i].kickValue;
		total += mult;
	}
	
	total /= (double)as.size();
	return total;
}

mat3x3 bond_angle_basis(mat3x3 parent, mat3x3 child)
{
	vec3 childDir = mat3x3_axis(child, 2);
	vec3 parentDir = mat3x3_axis(parent, 2);

}

double SVDBond::compareTorsions(BondPtr a, BondPtr b, bool ignore_intn)
{
	mat3x3 aBasis, bBasis;
	vec3 aPos, bPos;

	a->getAverageBasisPos(&aBasis, &aPos);
	b->getAverageBasisPos(&bBasis, &bPos);
	
	CorrelData cd = empty_CD();
	
	for (size_t i = 0; i < _atoms.size(); i++)
	{
		/* Exclude self, nan-potential */
		if (a->getMinor() == _atoms[i] || b->getMinor() == _atoms[i])
		{
			continue;
		}
		
		vec3 pos = _atoms[i]->getAbsolutePosition();
		
		/* The real important directions are rotated 90° but it doesn't
		 * really matter because we're comparing them */
		vec3 aDir = bond_effect_on_pos(pos, aBasis, aPos);
		vec3 bDir = bond_effect_on_pos(pos, bBasis, bPos);
		
		if (!ignore_intn)
		{
			double dir = _interactions[_atoms[i]][a];
			vec3_mult(&aDir, dir);

			dir = _interactions[_atoms[i]][b];
			vec3_mult(&bDir, dir);
		}
		
		for (int j = 0; j < 3; j++)
		{
			double *aVal = &aDir.x + j;
			double *bVal = &bDir.x + j;
			
			add_to_CD(&cd, *aVal, *bVal);
		}
	}
	
	double total = evaluate_CD(cd);

	return total;
}

double SVDBond::compareForKicks(BondPtr a, BondPtr b, bool kickmod)
{
	/* Get all the bond directions and positions */
	mat3x3 aBasis, bBasis;
	vec3 aPos, bPos;

	a->getAverageBasisPos(&aBasis, &aPos);
	b->getAverageBasisPos(&bBasis, &bPos);
	
	double total = 0;
	double count = 0;
	
	for (int i = 0; i < _atoms.size(); i++)
	{
		/* Exclude self, nan-potential */
		if (a->getMinor() == _atoms[i] || b->getMinor() == _atoms[i])
		{
			continue;
		}
		
		vec3 pos = _atoms[i]->getAbsolutePosition();
		
		/* The real important directions are rotated 90° but it doesn't
		 * really matter because we're comparing them */
		vec3 aDir = bond_effect_on_pos(pos, aBasis, aPos);
		vec3 bDir = bond_effect_on_pos(pos, bBasis, bPos);
		
		double aLength = vec3_length(aDir);
		double bLength = vec3_length(bDir);
		
		/* Get ratio of the lengths as one component of agreement */
		double ratio = aLength / bLength;
		if (ratio > 1)
		{
			ratio = 1 / ratio;
		}
		
		/* Get cosine of angle as another component of agreement */
		double cosine = vec3_cosine_with_vec3(aDir, bDir);
		
		double add = cosine * ratio;
		
		if (add != add)
		{
			continue;
		}
		
		total += add;
		
		count++;
	}
	
	total /= count;
	
	if (!kickmod)
	{
		return total;
	}

	double kicks = compareKicks(a, b);
	total *= kicks;

	return total;
}

void SVDBond::compareBonds()
{
	if (_doTorsion)
	{
		determineInteractions();
	}

	for (int i = 1; i < _bonds.size(); i++)
	{
		BondPtr b1 = _bonds[i];
		
		for (int j = 0; j < i; j++)
		{
			BondPtr b2 = _bonds[j];
			double agreement = 0;
			if (_doTorsion)
			{
				agreement = compareForKicks(b1, b2, false);
			}
			else
			{
				agreement = compareForKicks(b1, b2);
			}

			_svd[i][j] = agreement;
			_svd[j][i] = agreement;
		}
	}
	
	for (int i = 0; i < _bonds.size(); i++)
	{
		_svd[i][i] = 0;
	}
}

void SVDBond::performSVD()
{
	prepareMatrix(&_svd);
	prepareMatrix(&_original);
	prepareMatrix(&_v);
	prepareVector(&_w);
	
	compareBonds();
	
	writeMatrix();
	
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
	
	double pc50 = cuml * 0.45;
	
	if (_doTorsion)
	{
		pc50 = cuml * 0.9;
	}

	cuml = 0;
	_num = 0;
	for (int i = 0; i < _bonds.size(); i++)
	{
		cuml += _w[i];
		
		if (cuml > pc50)
		{
			_num = i;
			break;
		}
	}
	
	const int min = 10;
	if (_num < min)
	{
		_num = min;
		
		if (_bonds.size() < min)
		{
			_num = _bonds.size();
		}
	}
	
	/* Set up parameters for the top clusters from SVD */
	for (int i = 0; i < _num; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			ParamSVDSet set;
			Param *pw = NULL;
			Param *pk = NULL;

			if (j == 0)
			{
				pw = new Param();
				Param::setValue(pw, 0);
			}

			if (j == 1)
			{
				pk = new Param();
				Param::setValue(pk, 0);
			}

			Param *pph = new Param();
			Param::setValue(pph, 0);

			Param *t = new Param();
			Param::setValue(t, 0);

			set.pWhack = pw;
			set.pKick = pk;
			set.pPhi = pph;
			set.pTorsion = t;
			set.rowPtr = _svd[i];

			_params.push_back(set);
			
			if (j == 0 && _doTorsion)
			{
				break;
			}
		}
	}
	
	/* Get the baseline values for kicks and whacks */

	for (int i = 0; i < _bonds.size(); i++)
	{
		BondPtr bond = _bonds[i];
		double w = 0;
		double k = 0;

		if (bond->hasWhack())
		{
			WhackPtr whack = bond->getWhack();
			k = Whack::getKick(&*whack);
			w = Whack::getWhack(&*whack);
		}
		
		double t = Bond::getTorsion(&*bond);
		
		BondParamPair pair;
		pair.whack = w;
		pair.kick = k;
		pair.phi = 0;
		pair.torsion = t;
		pair.tr_x = t;
		pair.tr_y = t;
		pair.tr_z = t;
		_bondBase[bond] = pair;
	}
}

void SVDBond::addToStrategy(RefinementStrategyPtr strategy, double mult,
                            bool phi)
{
	double inv = 1.0 / _wTotal;
	inv = deg2rad(1.0);
	double tol = inv / 100;
	
	inv *= mult;
	
	for (int i = 0; i < _params.size() && !_doTorsion; i++)
	{
		strategy->addParameter(_params[i].pWhack, Param::getValue,
		                       Param::setValue, inv, tol);
		strategy->addParameter(_params[i].pKick, Param::getValue,
		                       Param::setValue, inv, tol);
	}
	
	for (int i = 0; i < _params.size() && _doTorsion; i++)
	{
		strategy->addParameter(_params[i].pTorsion, Param::getValue,
		                       Param::setValue, 
		                       deg2rad(mult), deg2rad(mult / 20));
	}
	
	double phi_step = deg2rad(180.);
	
	for (int i = 0; i < _params.size() && phi; i++)
	{
		strategy->addParameter(_params[i].pPhi, Param::getValue,
		                       Param::setValue, phi_step, deg2rad(5.));
	}
}

void SVDBond::bondsFromStrategy(RefinementStrategyPtr strategy)
{
	for (int i = 0; i < strategy->parameterCount(); i++)
	{
		Parameter obj = strategy->getParamObject(i);
		
		if (obj.getter == Bond::getTorsion)
		{
			BondPtr b = ((Bond *)(obj.object))->shared_from_this();
			_bonds.push_back(b);
		}
		
		if (obj.getter == Anchor::getPosX)
		{
			ModelPtr m = ((Model *)(obj.object))->shared_from_this();
			AnchorPtr a = ToAnchorPtr(m);
			_anchors.push_back(a);
		}
	}
}

void SVDBond::convertStrategyTorsions(RefinementStrategyPtr strategy,
                                      double torsion)
{
	for (int i = 0; i < strategy->parameterCount(); i++)
	{
		Parameter obj = strategy->getParamObject(i);
		
		if (obj.getter == Bond::getTorsion ||
		obj.getter == Anchor::getPosX ||
		obj.getter == Anchor::getPosY ||
		obj.getter == Anchor::getPosZ)
		{
			strategy->removeParameter(i);
			i--;
		}
	}

	for (int i = 0; i < _params.size(); i++)
	{
		strategy->addParameter(_params[i].pTorsion, Param::getValue,
		                       Param::setValue, 
		                       deg2rad(torsion), deg2rad(torsion / 20));
	}
}

void SVDBond::applyParameters()
{
	if (!_bonds.size())
	{
		return;
	}

	PolymerPtr pol = ToPolymerPtr(_bonds[0]->getAtom()->getMolecule());

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
		double ph = _bondBase[bi].phi;
		double t = _bondBase[bi].torsion;
		
		/* to store the total contributions to add to whack/kick */
		double w_more = 0;
		double k_more = 0;
		double ph_more = 0;
		double t_more = 0;

		/* loop round each cluster... */
		for (int j = 0; j < _params.size(); j++)
		{
			/* get the relative portion of this bond's identification
			 * with this cluster */
			
			double total = _params[j].rowPtr[i];
			
			/* get value of the whack/kick */
			if (!_doTorsion)
			{
				double wVal = Param::getValue(_params[j].pWhack);
				double kVal = Param::getValue(_params[j].pKick);
				double phVal = Param::getValue(_params[j].pPhi);

				wVal *= total;
				kVal *= total;
				phVal *= total;

				w_more += wVal;
				k_more += kVal;
				ph_more += phVal;
			}
			else
			{
				double tVal = Param::getValue(_params[j].pTorsion);

				if (!_silent)
				{
					std::cout << tVal << " " << total << std::endl;
				}

				tVal *= total;
				t_more += tVal;
			}

			if (k_more != k_more)
			{
				k_more = 0;
			}
			
			if (w_more != w_more)
			{
				w_more = 0;
			}
			
			if (ph_more != ph_more)
			{
				ph_more = 0;
			}
			
			if (t_more != t_more)
			{
				t_more = 0;
			}
		}

		if (w != w)
		{
			w = 0;
		}
		
		w += w_more;
		k += k_more;
		ph += ph_more;
		t += t_more;

		if (!_doTorsion)
		{
			WhackPtr whack = bi->getWhack();
			Whack::setWhack(&*whack, w);
			Whack::setKick(&*whack, k);
		}
		else
		{
			Bond::setTorsion(&*bi, t);
		}
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
		delete _params[i].pPhi;
		delete _params[i].pKick;
		delete _params[i].pTorsion;
	}
}

void SVDBond::writeMatrix()
{
	if (_silent)
	{
		return;
	}
	CSV csv(3, "bi", "bj", "cc");

	for (int i = 0; i < _bonds.size(); i++)
	{
		for (int j = 0; j < _bonds.size(); j++)
		{
			csv.addEntry(3, (double)i, (double)j, _svd[i][j]);
		}
	}

	csv.setSubDirectory("local_flex");
	csv.writeToFile("new_bond_matrix.csv");

	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = "new_bond_matrix";
	plotMap["height"] = "1000";
	plotMap["width"] = "1000";
	plotMap["xHeader0"] = "bi";
	plotMap["yHeader0"] = "bj";
	plotMap["zHeader0"] = "cc";

	plotMap["xTitle0"] = "bond number";
	plotMap["yTitle0"] = "bond number";

	plotMap["style0"] = "heatmap";
	plotMap["stride"] = i_to_str(_bonds.size());
	
	csv.setSubDirectory("local_flex");
	csv.plotPNG(plotMap);
}

void SVDBond::interactionsForBond(BondPtr bond)
{
	AtomGroupPtr downstream = bond->makeAtomGroup();
	
	for (size_t i = 0; i < downstream->atomCount(); i++)
	{
		AtomPtr a = downstream->atom(i);
		if (a->getElectronCount() <= 1)
		{
			continue;
		}
		
		if (std::find(_atoms.begin(), _atoms.end(), a) == _atoms.end())
		{
			if (!_silent)
			{
				std::cout << "Don't have " << a->shortDesc() << std::endl;
			}

			continue;
		}

		if (!_silent && false)
		{
			std::cout << bond->shortDesc() << " interacts with " << a->shortDesc() << std::endl;
		}

		_interactions[a][bond] = 1;
	}
	
	if (!bond->hasTwist())
	{
		return;
	}
	
	for (size_t i = 0; i < _atoms.size(); i++)
	{
		if (downstream->hasAtom(_atoms[i]))
		{
			continue;
		}

		AtomPtr a = _atoms[i];

		if (a->getElectronCount() <= 1)
		{
			continue;
		}
		
		if (bond->getMajor() == a || bond->getMinor() == a)
		{
			continue;
		}
		
		_interactions[a][bond] = -1;
	}
}

void SVDBond::determineInteractions()
{
	for (size_t i = 0; i < _bonds.size(); i++)
	{
		interactionsForBond(_bonds[i]);
	}
}
