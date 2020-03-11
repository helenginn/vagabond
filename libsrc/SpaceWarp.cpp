// Vagabond
// Copyright (C) 2019 Helen Ginn
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

#include "../libica/svdcmp.h"
#include "SpaceWarp.h"
#include "SVDBond.h"
#include "RefinementNelderMead.h"
#include "FFT.h"
#include "Timer.h"
#include "Atom.h"
#include "AtomGroup.h"
#include "Options.h"
#include "Bond.h"
#include "maths.h"
#include <map>

SpaceWarp::SpaceWarp(VagFFTPtr fft)
{
	_fft = fft;
	long nn = _fft->nn();
	
	_atomCounter = 0;
	_bondCounter = 0;
	
	_minRegion = empty_vec3();
	_maxRegion = empty_vec3();

	_warps = (vec3 *)malloc(nn * sizeof(vec3));
	memset(_warps, '\0', nn * sizeof(vec3));
}

vec3 SpaceWarp::target(AtomPtr atom)
{
	mat3x3 rb = _fft->toRecip();
	vec3 pos = atom->getPositionInUnitCell();
	mat3x3_mult_vec(rb, &pos);
	_fft->expandToVoxel(&pos);
	vec3 buffer = make_vec3(0, 0, 0);

	_minRegion = vec3_subtract_vec3(pos, buffer);
	_maxRegion = vec3_add_vec3(pos, buffer);
	
	long ele = _fft->element(_minRegion.x, _minRegion.y, _minRegion.z);

	return _warps[ele];
}

void SpaceWarp::initAtom(AtomPtr atom)
{
	mat3x3 rb = _fft->toRecip();
	vec3 pos = atom->getPositionInUnitCell();
	mat3x3_mult_vec(rb, &pos);
	_fft->expandToVoxel(&pos);
	vec3 buffer = make_vec3(2, 2, 2);

	_minRegion = vec3_subtract_vec3(pos, buffer);
	_maxRegion = vec3_add_vec3(pos, buffer);
	
	for (int z = _minRegion.z + 1; z <= _maxRegion.z - 1; z++)
	{
		for (int y = _minRegion.y + 1; y <= _maxRegion.y - 1; y++)
		{
			for (int x = _minRegion.x + 1; x <= _maxRegion.x - 1; x++)
			{
				long ele = _fft->element(x, y, z);
				
				std::vector<AnyPtr> anys;
				for (int i = 0; i < 3; i++)
				{
					AnyPtr any = AnyPtr(new Any(&_warps[ele].x + i));
					_varying.push_back(any);
				}
			}
		}
	}
}

double SpaceWarp::evaluate()
{
	double sum_x = 0;
	double sum_y = 0;
	double sum_xx = 0;
	double sum_yy = 0;
	double sum_xy = 0;
	double sum_w = 0;
	double penalty = 0;

	for (int z = _minRegion.z; z < _maxRegion.z; z++)
	{
		for (int y = _minRegion.y; y <= _maxRegion.y; y++)
		{
			for (int x = _minRegion.x; x <= _maxRegion.x; x++)
			{
				long ele = _fft->element(x, y, z);
				vec3 pos = make_vec3(x, y, z);
				vec3 nudge = _warps[ele];

				for (int l = -1; l <= 1; l++)
				{
					for (int k = -1; k <= 1; k++)
					{
						for (int h = -1; h <= 1; h++)
						{
							long mele = _fft->element(x + h, y + k, z + l);
							double mult = exp(-(h * h + k * k + l * l));
							vec3 mini_nudge = _warps[mele];
							vec3_mult(&mini_nudge, mult);

							vec3_add_to_vec3(&pos, nudge);
						}
					}
				}

				double nx = _data->cubic_interpolate(pos, 0);
				double y = _fft->getReal(ele);
				double weight = 1;

				sum_x += nx * weight;
				sum_y += y * weight;
				sum_yy += y * y * weight;
				sum_xx += nx * nx * weight;
				sum_xy += nx * y * weight;
				sum_w += weight;
			}
		}
	}

	double top = sum_w * sum_xy - sum_x * sum_y;
	double bottom_left = sum_w * sum_xx - sum_x * sum_x;
	double bottom_right = sum_w * sum_yy - sum_y * sum_y;
	
	double r = top / sqrt(bottom_left * bottom_right);
	
//	std::cout << correl << std::endl;
	return -r;
}

void SpaceWarp::nudgeAtom(AtomPtr atom)
{
	initAtom(atom);
	
	NelderMeadPtr neld = NelderMeadPtr(new RefinementNelderMead());
	neld->setVerbose(false);
	neld->setSilent(true);
	neld->setCycles(10);
	neld->setEvaluationFunction(SpaceWarp::getScore, this);
	
	for (int i = 0; i < _varying.size(); i += 3)
	{
		for (int j = 0; j < 3; j++)
		{
			AnyPtr any = _varying[i + j];
			neld->addParameter(&*any, Any::get, Any::set, 0.1, 0.0001,
			                   "a" + i_to_str(j));
		}
		neld->refine();
		neld->clearParameters();
	}
	

	_varying.clear();
}

void SpaceWarp::recalculate(VagFFTPtr data)
{
	_bonds.clear();
	_whacks.clear();
	_atoms = AtomGroupPtr();
	_atomBonds.clear();

	Options::getRuntimeOptions()->renderDensity();
	long nn = _fft->nn();
	memset(_warps, '\0', nn * sizeof(vec3));

	Timer timer;
	std::cout << "Calculating local map warp... " << std::flush;
	_data = data;
	_fft->fft(FFTReciprocalToReal);
	mat3x3 rb = _fft->getRealBasis();
	
	CrystalPtr crystal = Options::getActiveCrystal();
	
	for (int i = 0; i < crystal->atomCount(); i++)
	{
		AtomPtr a = crystal->atom(i);
		if (a->isBackbone() && a->getElectronCount() > 1)
		{
			nudgeAtom(a);
		}
	}

	Options::getRuntimeOptions()->renderWarp();

	_fft->fft(FFTRealToReciprocal);

	for (int i = 0; i < _fft->nn(); i++)
	{
		mat3x3_mult_vec(rb, &_warps[i]);
	}
	
	timer.quickReport();
	std::cout << std::endl;
}

SpaceWarp::~SpaceWarp()
{
	long nn = _fft->nn();
	memset(_warps, '\0', nn * sizeof(vec3));
	free(_warps);
	_warps = NULL;
}

void SpaceWarp::addRefinedAtom(AtomPtr atom)
{
	if (!atom->getModel()->isBond())
	{
		return;
	}

	BondPtr bond = ToBondPtr(atom->getModel());
	addBondAtoms(bond);

	if (!bond->hasWhack())
	{
		return;
	}
	
	addWhack(bond->getWhack());
}

void SpaceWarp::addBondAtoms(BondPtr bond)
{
	addBond(bond);
	
	if (!_atoms)
	{
		_atoms = AtomGroupPtr(new AtomGroup());
	}

	AtomGroupPtr atoms = bond->makeAtomGroup();
	
	for (int i = 0; i < atoms->atomCount(); i++)
	{
		AtomPtr a = atoms->atom(i);
		if (!a->isBackbone() || a->getElectronCount() <= 1)
		{
			continue;
		}

		if (!_atomBonds.count(a))
		{
			_atomBonds[a] = BondList();
		}

		_atomBonds[a].push_back(bond);
		_atoms->addAtom(a);
	}
}

vec3 SpaceWarp::bondEffect(vec3 atomPos, BondPtr b)
{
	mat3x3 bBasis;
	vec3 bPos;
	b->getAverageBasisPos(&bBasis, &bPos);
	vec3 derivative = bond_effect_on_pos(atomPos, bBasis, bPos);

	return derivative;
}

void SpaceWarp::populateSVD(BondPtr bond)
{
	mat3x3 toRealPos = _fft->getRealBasis();

	vec3 xyz = make_vec3(_minRegion.x, _minRegion.y, _minRegion.z);
	mat3x3_mult_vec(toRealPos, &xyz);
	vec3 effect = bondEffect(xyz, bond);

	_matPtrs[_bondCounter][_atomCounter + 0] = effect.x;
	_matPtrs[_bondCounter][_atomCounter + 1] = effect.y;
	_matPtrs[_bondCounter][_atomCounter + 2] = effect.z;
}

void SpaceWarp::populateSVD(AtomPtr atom, BondList bonds)
{
	mat3x3 rb = _fft->toRecip();
	vec3 pos = atom->getPositionInUnitCell();
	mat3x3_mult_vec(rb, &pos);
	_fft->expandToVoxel(&pos);
	vec3 buffer = make_vec3(0, 0, 0);
	
	_minRegion = vec3_subtract_vec3(pos, buffer);
	_maxRegion = vec3_add_vec3(pos, buffer);
	
	for (int i = 0; i < bonds.size(); i++)
	{
		BondPtr b = bonds[i];
		populateSVD(b);
		_bondCounter++;
	}
}

void SpaceWarp::addTargets()
{
	for (int i = 0; i < _atoms->atomCount(); i++)
	{
		vec3 t = target(_atoms->atom(i));
		
		for (int j = 0; j < 3; j++)
		{
			int index = i * 3 + j;
			_t[index] = *(&t.x + j);
		}
	}
}

void SpaceWarp::svd()
{
	setupSVD();
	
	_atomCounter = 0;

	for (int i = 0; i < _atoms->atomCount(); i++)
	{
		AtomPtr a = _atoms->atom(i);
		BondList bonds = _atomBonds[a];
		
		populateSVD(a, bonds);
		_atomCounter += 3;
		_bondCounter = 0;
	}
	
	addTargets();

	int bDim = _bonds.size();
	int aDim = _atoms->atomCount() * 3;

	int success = svdcmp((mat)_matPtrs, aDim, bDim, (vect)_w, (mat)_vPtrs);
	
	std::cout << "Success: " << success << std::endl;
	
	if (!success)
	{
		std::cout << "shit." << std::endl;
		cleanupSVD();
		return;
	}
	
	double *tmp = (double *)malloc(sizeof(double *) * bDim);

	double sum = 0;
	for (int j = 0; j < bDim; j++)
	{
		sum += _w[j]; 
	}
	
	sum /= bDim;

	for (int j = 0; j < bDim; j++)
	{
		if (_w[j] < sum * 0.1)
		{
			continue;
		}

		double sum = 0;
		for (int i = 0; i < aDim; i++)
		{
			sum += _matPtrs[j][i] * _t[i];
		}
		
		tmp[j] = 1 / _w[j] * sum;
	}
	
	double step = 0.2;
	
	std::cout << "Weights: " << std::endl;
	for (int i = 0; i < bDim; i++)
	{
		double sum = 0;
		for (int j = 0; j < bDim; j++)
		{
			sum += _vPtrs[i][j] * tmp[j];
		}
		
		_weights[i] = sum * step;
		std::cout << _weights[i] << std::endl;
	}
	
	for (int i = 0; i < bDim; i++)
	{
		BondPtr b = _bonds[i];
		double t = Bond::getTorsion(&*b);
		t += _weights[i];
		Bond::setTorsion(&*b, t);
	}
	
	Options::getActiveCrystal()->refreshPositions();
	
	free(tmp);

	cleanupSVD();
}

void SpaceWarp::setupSVD()
{
	size_t atomCount = _atoms->atomCount() * 3; /* 3 = 3D */
	size_t bondCount = _bonds.size();
	double rectsize = sizeof(double) * atomCount * bondCount;
	double squaresize = sizeof(double) * bondCount * bondCount;

	_matrix = (double *)malloc(rectsize);
	_v = (double *)malloc(squaresize);

	memset(_matrix, '\0', rectsize);
	memset(_v, '\0', squaresize);
	
	_matPtrs = (double **)malloc(sizeof(double **) * atomCount);
	_vPtrs = (double **)malloc(sizeof(double **) * bondCount);
	_w = (double *)malloc(sizeof(double *) * bondCount);
	_t = (double *)malloc(sizeof(double *) * atomCount);
	_weights = (double *)malloc(sizeof(double *) * bondCount);

	memset(_w, '\0', sizeof(double) * bondCount);
	memset(_t, '\0', sizeof(double) * atomCount);
	memset(_weights, '\0', sizeof(double) * bondCount);
	
	for (int i = 0; i < atomCount; i++)
	{
		_matPtrs[i] = &_matrix[i * bondCount];
	}

	for (int i = 0; i < bondCount; i++)
	{
		_vPtrs[i] = &_v[i * bondCount];
	}
}

void SpaceWarp::cleanupSVD()
{
	free(_matrix);
	free(_w);
	free(_v);
	free(_matPtrs);
	free(_vPtrs);
	free(_weights);
	
	_matrix = NULL;
	_w = NULL;
	_v = NULL;
	_matPtrs = NULL;
	_vPtrs = NULL;
	_weights = NULL;
}
