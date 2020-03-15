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
#include "CSV.h"
#include "Timer.h"
#include "Twist.h"
#include "Atom.h"
#include "AtomGroup.h"
#include "FlexGlobal.h"
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
	
	_activePos = empty_vec3();

	_warps = (vec3 *)malloc(nn * sizeof(vec3));
	memset(_warps, '\0', nn * sizeof(vec3));
}

vec3 SpaceWarp::target(AtomPtr atom)
{
	mat3x3 rb = _fft->toRecip();
	vec3 pos = atom->getPositionInUnitCell();
	mat3x3_mult_vec(rb, &pos);
	_fft->expandToVoxel(&pos);

	_activePos = pos;
	
	long ele = _fft->element(_activePos.x, _activePos.y, _activePos.z);

	return _warps[ele];
}

void SpaceWarp::initAtom(AtomPtr atom)
{
	mat3x3 rb = _fft->toRecip();
	vec3 pos = atom->getPositionInUnitCell();
	mat3x3_mult_vec(rb, &pos);
	_fft->expandToVoxel(&pos);

	_activePos = pos;
	
	long ele = _fft->element(pos.x, pos.y, pos.z);
	_warps[ele] = empty_vec3();

	std::vector<AnyPtr> anys;
	for (int i = 0; i < 3; i++)
	{
		AnyPtr any = AnyPtr(new Any(&_warps[ele].x + i));
		_varying.push_back(any);
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
	
	vec3 c = _activePos;

	long ele = _fft->element(c.x, c.y, c.z);
	vec3 nudge = _warps[ele];
	int buf = 3;

	for (int l = -buf; l <= buf; l++)
	{
		for (int k = -buf; k <= buf; k++)
		{
			for (int h = -buf; h <= buf; h++)
			{
				vec3 nbr = make_vec3(c.x + h, c.y + k, c.z + l);
				vec3 shift = vec3_add_vec3(nbr, nudge);

				double x = _fft->cubic_interpolate(nbr, 0);
				double y = _data->cubic_interpolate(shift, 0);
				double weight = 1;

				sum_x += x * weight;
				sum_y += y * weight;
				sum_yy += y * y * weight;
				sum_xx += x * x * weight;
				sum_xy += x * y * weight;
				sum_w += weight;
			}
		}
	}
	
	mat3x3 toRealPos = _fft->getRealBasis();
	mat3x3_mult_vec(toRealPos, &nudge);
	double penalty = exp(-vec3_sqlength(nudge));

	double top = sum_w * sum_xy - sum_x * sum_y;
	double bottom_left = sum_w * sum_xx - sum_x * sum_x;
	double bottom_right = sum_w * sum_yy - sum_y * sum_y;
	
	double r = top / sqrt(bottom_left * bottom_right);
	r *= penalty;
	
	return -r;
}

void SpaceWarp::nudgeAtom(AtomPtr atom)
{
	initAtom(atom);
	
	NelderMeadPtr neld = NelderMeadPtr(new RefinementNelderMead());
	neld->setVerbose(false);
	neld->setSilent(true);
	neld->setCycles(200);
	neld->setEvaluationFunction(SpaceWarp::getScore, this);
	mat3x3 rb = _fft->getRealBasis();
	
	vec3 sum = empty_vec3();
	vec3 pos = _activePos;
	long ele = _fft->element(pos.x, pos.y, pos.z);

	for (int j = 0; j < 3; j++)
	{
		double step = 0.5;

		AnyPtr any = _varying[j];
		neld->addParameter(&*any, Any::get, Any::set, step, 0.0001,
		                   atom->shortDesc() + "_" + i_to_str(j));
	}

	neld->refine();

	vec3 warp = _warps[ele];

	neld->clearParameters();

	_varying.clear();
}

void SpaceWarp::nudge()
{
	long nn = _fft->nn();
	_fft->wipe();
	_atoms->addToMap(_fft);

	mat3x3 rb = _fft->getRealBasis();
	
	for (int i = 0; i < _atoms->atomCount(); i++)
	{
		AtomPtr a = _atoms->atom(i);
		
		if (a->getElectronCount() > 1)
		{
			nudgeAtom(a);
		}
	}

	Options::getRuntimeOptions()->renderWarp();

	for (int i = 0; i < _fft->nn(); i++)
	{
		mat3x3_mult_vec(rb, &_warps[i]);
	}

}

void SpaceWarp::recalculate(VagFFTPtr data)
{
	Timer timer;
	std::cout << "Nudging atoms around... " << std::flush;

	_bonds.clear();
	_interactions.clear();
	_data = data;

	for (int i = 0; i < _atoms->atomCount(); i++)
	{
		AtomPtr a = _atoms->atom(i);
		if (a->getModel()->isBond() && 
		    ToBondPtr(a->getModel())->isRefinable() && 
		    ToBondPtr(a->getModel())->isUsingTorsion())
		{

			addRefinedAtom(a);
		}
	}
	

	setupSVD();
	for (int i = 0; i < 1; i++)
	{
		nudge();
		svd();
	}

	cleanupSVD();
	
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
}

void SpaceWarp::addBondAtoms(BondPtr bond)
{
	addBond(bond);
	
	AtomGroupPtr downstream = bond->makeAtomGroup();
	
	for (int i = 0; i < downstream->atomCount(); i++)
	{
		AtomPtr a = downstream->atom(i);
		if (a->getElectronCount() <= 1)
		{
			continue;
		}

		_interactions[a][bond] = 1;

		_atomsForSVD->addAtom(a);
	}
	
	if (!bond->hasTwist())
	{
		return;
	}
	
	for (int i = 0; i < _atoms->atomCount(); i++)
	{
		if (downstream->hasAtom(_atoms->atom(i)))
		{
			continue;
		}

		AtomPtr a = _atoms->atom(i);

		if (a->getElectronCount() <= 1)
		{
			continue;
		}
		
		_interactions[a][bond] = -1;

		_atomsForSVD->addAtom(a);
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

void SpaceWarp::populateSVD(AtomPtr a, BondPtr bond, int interaction)
{
	vec3 atomPos = a->getAbsolutePosition();
	vec3 effect = bondEffect(atomPos, bond);

	_matPtrs[_atomCounter + 0][_bondCounter] = effect.x * abs(interaction);
	_matPtrs[_atomCounter + 1][_bondCounter] = effect.y * abs(interaction);
	_matPtrs[_atomCounter + 2][_bondCounter] = effect.z * abs(interaction);
}

void SpaceWarp::populateSVD(AtomPtr atom)
{
	mat3x3 rb = _fft->toRecip();
	vec3 pos = atom->getPositionInUnitCell();
	mat3x3_mult_vec(rb, &pos);
	_fft->expandToVoxel(&pos);
	vec3 buffer = make_vec3(0, 0, 0);
	
	_activePos = pos;
	
	_bondCounter = 0;

	for (int i = 0; i < _bonds.size(); i++)
	{
		BondPtr b = _bonds[i];
		int interaction = 0;
		
		if (_interactions.count(atom) && _interactions[atom].count(b))
		{
			interaction = _interactions[atom][b];
			
			if (interaction > 0)
			{
				populateSVD(atom, b, interaction);
			}
		}

		_bondCounter++;
		
		if (_interactions.count(atom) && _interactions[atom].count(b))
		{
			interaction = _interactions[atom][b];
			
			if (interaction < 0)
			{
				populateSVD(atom, b, interaction);
			}
		}

		_bondCounter++;
	}
}

void SpaceWarp::addTargets()
{
	double sum = 0;
	for (int i = 0; i < _atomsForSVD->atomCount(); i++)
	{
		vec3 t = target(_atomsForSVD->atom(i));
		
		for (int j = 0; j < 3; j++)
		{
			int index = i * 3 + j;
			_t[index] = -*(&t.x + j);
			sum += *(&t.x + j) * *(&t.x + j);
		}
	}
	std::cout << "Target sum: " << sum << std::endl;
}

void SpaceWarp::svd()
{
	std::cout << "Starting SVD..." << std::endl;
	Timer t;
	
	_atomCounter = 0;

	for (int i = 0; i < _atomsForSVD->atomCount(); i++)
	{
		AtomPtr a = _atomsForSVD->atom(i);
		populateSVD(a);

		_atomCounter += 3;
	}
	
	addTargets();

	int bDim = _bonds.size() * 2;
	int aDim = _atomsForSVD->atomCount() * 3;
	
	CSVPtr csv = CSVPtr(new CSV(3, "b", "a", "v"));
	csv->setSubDirectory("local_flex");
	
	for (int j = 0; j < bDim; j++)
	{
		for (int i = 0; i < aDim; i++)
		{
			csv->addEntry(3, (double)j, (double)i, (double)(_matPtrs[i][j]));
		}
	}

	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = "spacewarp";
	plotMap["height"] = "3400";
	plotMap["width"] = "1400";
	plotMap["xHeader0"] = "b";
	plotMap["yHeader0"] = "a";
	plotMap["zHeader0"] = "v";

	plotMap["xTitle0"] = "b dim";
	plotMap["yTitle0"] = "a dim";
	plotMap["style0"] = "heatmap";
	plotMap["stride"] = i_to_str(bDim);

	csv->writeToFile("spacewarp.csv");
	csv->plotPNG(plotMap);

	csv = CSVPtr(new CSV(3, "b", "a", "v"));
	csv->setSubDirectory("local_flex");
	
	for (int j = 0; j < _bonds.size(); j++)
	{
		for (int i = 0; i < _atomsForSVD->atomCount(); i++)
		{
			double val = _interactions[_atomsForSVD->atom(i)][_bonds[j]];
			csv->addEntry(3, (double)j, (double)i, val);
		}
	}

	plotMap["stride"] = i_to_str(_bonds.size());
	plotMap["filename"] = "interactions";
	csv->plotPNG(plotMap);


	std::cout << "Doing the SVD..." << std::endl;
	int success = svdcmp((mat)_matPtrs, aDim, bDim, (vect)_w, (mat)_vPtrs);
	
	if (!success)
	{
		std::cout << "shit." << std::endl;
		std::cout << "m x n = " << aDim << " " << bDim << std::endl;
		return;
	}
	
	double *tmp = (double *)malloc(sizeof(double *) * bDim);
	memset(tmp, '\0', sizeof(double *) * bDim);

	double sum = 0;
	for (int j = 0; j < bDim; j++)
	{
		sum += _w[j]; 
	}
	
	sum /= bDim;

	for (int j = 0; j < bDim; j++)
	{
		BondPtr b = _bonds[j / 2];
		if (_w[j] < sum * 0.001)
		{
//			if (b->makeAtomGroup()->atomCount() > 1)
			{
				continue;
			}
		}

		double sum = 0;
		for (int i = 0; i < aDim; i++)
		{
			sum += _matPtrs[i][j] * _t[i];
		}
		
		tmp[j] = 1 / _w[j] * sum;
		
		if (tmp[j] != tmp[j] || !std::isfinite(tmp[j]))
		{
			tmp[j] = 0;
		}
	}
	
	for (int i = 0; i < bDim; i++)
	{
		double sum = 0;
		for (int j = 0; j < bDim; j++)
		{
			sum += _vPtrs[i][j] * tmp[j];
		}
		
		_weights[i] = sum;
	}
	
	for (int i = 0; i < _bonds.size(); i++)
	{
		BondPtr b = _bonds[i];
		double t = Bond::getTorsion(&*b);
		_torsions[i * 2] = t;
	}

	_magnitude = 10 / (double)_bonds.size();
	
	NelderMeadPtr grid;
	grid = NelderMeadPtr(new RefinementNelderMead());
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	_target = new FlexGlobal();
	_target->setCrystal(crystal);
	_target->setAtomGroup(_all);
	AnyPtr any = AnyPtr(new Any(&_magnitude));
	grid->addParameter(&*any, Any::get, Any::set, 
	                   _magnitude * 0.2, _magnitude * 0.0001);
	_magnitude = 0;
//	_magnitude = 0.02;
	grid->setVerbose(true);
	grid->setSilent(false);
	grid->setCycles(100);
	grid->setEvaluationFunction(SpaceWarp::staticScore, this);
	FlexGlobal::score(_target);
	score();
	grid->refine();
	
	delete _target;
	free(tmp);
	t.quickReport();
}

double SpaceWarp::score()
{
	int bDim = _bonds.size();
	for (int i = 0; i < bDim; i++)
	{
		BondPtr b = _bonds[i];
		double t = _torsions[i * 2];
		t += _weights[i * 2] * _magnitude;
		Bond::setTorsion(&*b, t);
		
		if (b->hasTwist())
		{
			double tw = _weights[i * 2 + 1] * _magnitude;
			Twist::setTwist(&*b->getTwist(), tw);
			b->getTwist()->getAppliedModel()->propagateChange();
		}
	}
	
	_atoms->refreshPositions(false);
	return FlexGlobal::score(_target);
}

void SpaceWarp::setupSVD()
{
	_atomsForSVD = _atoms;
	size_t atomCount = _atomsForSVD->atomCount() * 3;
	size_t bondCount = _bonds.size() * 2;
	double rectsize = sizeof(double) * atomCount * bondCount;
	double squaresize = sizeof(double) * bondCount * bondCount;

	_matrix = (double *)malloc(rectsize);
	_v = (double *)malloc(squaresize);
	
	_matPtrs = (double **)malloc(sizeof(double **) * atomCount);
	_vPtrs = (double **)malloc(sizeof(double **) * bondCount);
	_w = (double *)malloc(sizeof(double *) * bondCount);
	_t = (double *)malloc(sizeof(double *) * atomCount);
	_weights = (double *)malloc(sizeof(double *) * bondCount);

	_torsions = (double *)malloc(sizeof(double *) * bondCount);
	
	for (int i = 0; i < atomCount; i++)
	{
		_matPtrs[i] = &_matrix[i * bondCount];
	}

	for (int i = 0; i < bondCount; i++)
	{
		_vPtrs[i] = &_v[i * bondCount];
	}
	
	wipeSVD();
}

void SpaceWarp::wipeSVD()
{
	size_t atomCount = _atomsForSVD->atomCount() * 3;
	size_t bondCount = _bonds.size() * 2;

	memset(_matrix, '\0', sizeof(double) * atomCount * bondCount);
	memset(_v, '\0', sizeof(double) * bondCount * bondCount);
	memset(_w, '\0', sizeof(double) * bondCount);
	memset(_t, '\0', sizeof(double) * atomCount);
	memset(_weights, '\0', sizeof(double) * bondCount);
	memset(_torsions, '\0', sizeof(double) * bondCount);
}

void SpaceWarp::cleanupSVD()
{
	free(_matrix);
	free(_w);
	free(_v);
	free(_matPtrs);
	free(_vPtrs);
	free(_weights);
	free(_torsions);
	
	_matrix = NULL;
	_w = NULL;
	_v = NULL;
	_matPtrs = NULL;
	_vPtrs = NULL;
	_weights = NULL;
	_torsions = NULL;
}

void SpaceWarp::setCalculated(VagFFTPtr fft)
{
	_fft = VagFFTPtr(new VagFFT(*fft));
	_fft->fft(FFTReciprocalToReal);
}

void SpaceWarp::setActiveAtoms(AtomGroupPtr atoms)
{
	_atoms = AtomGroupPtr(new AtomGroup());
	_atomsForSVD = AtomGroupPtr(new AtomGroup());
	_all = atoms;
	
	for (int i = 0; i < atoms->atomCount(); i++)
	{
		AtomPtr a = atoms->atom(i);
		
		if (a->getModel()->isAnchor())
		{
			_atoms->addAtom(a);
		}

		if (a->getAtomName() != "CA")
		{
			continue;
		}

		if (a->getElectronCount() > 1)
		{
			_atoms->addAtom(a);
		}
	}
}
