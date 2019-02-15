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

#include <iomanip>
#include "Gradiator.h"
#include "Anchor.h"
#include "Crystal.h"
#include "Polymer.h"
#include "Options.h"
#include "Bond.h"
#include "Atom.h"
#include "Bucket.h"
#include "fftw3d.h"
#include "Whack.h"

Gradiator::Gradiator(PolymerPtr polymer)
{
	_crystal = Options::getActiveCrystal();
	_polymer = polymer;
}

bool in_asu(CSym::CCP4SPG *spg, vec3 vec)
{
	if (vec.x < spg->mapasu_zero[0] &&
	    vec.y < spg->mapasu_zero[1] &&
	    vec.z < spg->mapasu_zero[2]) 
	{
		return true;
	}
	
	return false;
}

void Gradiator::prepareList()
{
	AtomGroupPtr backbone = _polymer->getAllBackbone();
	FFTPtr fft = _crystal->getFFT();
	BucketPtr solv = _crystal->getBucket();

	double cutoff_dsq = 4.0;
	double all_atoms = 0;
	double plausible_voxels = 0;
	mat3x3 f2rt = _crystal->getReal2Frac();
	mat3x3 f2rt2 = mat3x3_inverse(f2rt);
	mat3x3 f2r = mat3x3_transpose(f2rt2);
	CSym::CCP4SPG *spg = _crystal->getSpaceGroup();
	
	for (int i = 0; i < backbone->atomCount(); i++)
	{
		AtomPtr a = backbone->atom(i);

		if (!a->getModel()->hasExplicitPositions())
		{
			continue;
		}
		
		vec3 asu = a->getPositionInAsu();
		_atoms[a] = asu;
	}
	
	for (size_t j = 0; j < fft->nn; j++)
	{
		vec3 real = fft->fracFromElement(j);
		
		if (!in_asu(spg, real))
		{
			continue;
		}
		
		plausible_voxels++;

		if (solv && solv->isSolvent(j))
		{
			continue;
		}

		mat3x3_mult_vec(f2r, &real);
		
		Voxel vox;
		vox.pos = real;
		vox.obs = fft->data[j][0];
		vox.calc = fft->data[j][1];

		std::vector<SingleAtom> nearby;
		
		for (int i = 0; i < backbone->atomCount(); i++)
		{
			AtomPtr a = backbone->atom(i);
			
			if (_atoms.count(a) == 0)
			{
				continue;
			}
			
			vec3 atom_real = _atoms[a];
			vec3 diff = vec3_subtract_vec3(atom_real, real);

			if (vec3_sqlength(diff) < cutoff_dsq)
			{
				ExplicitModelPtr m = a->getExplicitModel();
				std::vector<BondSample> *samples = m->getManyPositions();
				
				nearby.reserve(nearby.size() + samples->size());
				for (int k = 0; k < samples->size(); k++)
				{
					SingleAtom pair;
					pair.index = k;
					vec3 start = m->getManyPositions()->at(k).start;
					vec3 end = FFT::getPositionInAsu(start);
					pair.pos = end;
					pair.mag = m->getManyPositions()->at(k).kickValue;
					nearby.push_back(pair);
					all_atoms++;
				}
			}
		}
		
		if (nearby.size())
		{
			vox.nearby_atoms = nearby;
			_voxels.push_back(vox);
		}
	}
	
	std::cout << "Total: " << _voxels.size() << " out of " 
	<< plausible_voxels << std::endl;
	
	std::cout << "Atoms per voxel: " << all_atoms / (double)_voxels.size()
	<< std::endl;
	
	AnchorPtr anchor = _polymer->getAnchorModel();
	
	for (int i = 0; i < anchor->whackCount(); i++)
	{
		WhackPtr whack = anchor->getWhack(i);
		WhackVal wv;
		wv.change = 0;
		BondPtr bond = whack->getBond();
		
		std::vector<BondSample> *samples = bond->getManyPositions();

		for (int i = 0; i < samples->size(); i++)
		{
			wv.pos.push_back(samples->at(i).start);
			vec3 dir = mat3x3_axis(samples->at(i).basis, 2);
			wv.dir.push_back(dir);
		}

		_ws.push_back(whack);
		_whacks[whack] = wv;
	}
	
	std::cout << "Whacks: " << _ws.size() << std::endl;
}

double Gradiator::dDistanceTodDensity(AtomPtr a, double curr, double delta)
{
	curr *= 2;
	double grad = -2 * curr * exp(-curr * curr);
	grad *= delta;

	return grad;
}

double deltaSqrtFunc(double val)
{
	double ret = 0.5 / sqrt(val);
	
	if (ret != ret)
	{
		ret = 0;
	}
	
	return ret;
}

vec3 Gradiator::whackPosition(Voxel *vox, size_t n, WhackPtr whack)
{
	std::vector<SingleAtom> &nearby = vox->nearby_atoms;
	size_t index = nearby[n].index;
	
	if (_whacks.count(whack))
	{
		return _whacks[whack].pos[index];
	}
	
	return empty_vec3();
}

double currentDistance(Voxel *vox, size_t n)
{
	vec3 atom_pos = vox->nearby_atoms[n].pos;
	vec3_subtract_from_vec3(&atom_pos, vox->pos);
	
	double length = vec3_length(atom_pos);
	return length;
}

double Gradiator::deltaDistance(Voxel *vox, size_t n, 
                                WhackPtr whack)
{
	double mag = vox->nearby_atoms[n].mag;
	int index = vox->nearby_atoms[n].index;
	vec3 vp = vox->pos;
	vec3 bp = _whacks[whack].pos[index];
	vec3 bd = _whacks[whack].dir[index];
	double change = _whacks[whack].change;
	vec3 ap = vox->nearby_atoms[n].pos;
	vec3 sub = vec3_subtract_vec3(ap, bp);
	
	double xx = bd.y * sub.z - bd.z * sub.y;
	xx *= mag;
	double xa = ap.x - vp.x + change * xx;

	double yx = bd.z * sub.x - bd.x * sub.z;
	yx *= mag;
	double ya = ap.y - vp.y + change * yx;

	double zx = bd.x * sub.y - bd.y * sub.x;
	zx *= mag;
	double za = ap.z - vp.z + change * zx;

	/* chained derivative */
	double add = xa * xa + ya * ya + za * za;
	add = 0.5 / sqrt(add);

	double d_mult = 2 * xa * xx + 2 * ya * yx + 2 * za * zx;
	double deriv = add * d_mult;
	
	return deriv;
}

double Gradiator::deltaVoxelValue(Voxel *vox)
{
	double cuml = 0;

//	...
	
	return cuml;
}

double Gradiator::deltaVoxel4Whack(WhackPtr whack, Voxel *vox)
{
	std::vector<SingleAtom> *nearby = &vox->nearby_atoms;
	
	double cuml = 0;
	
	for (size_t i = 0; i < nearby->size(); i++)
	{
		double delta_dist = deltaDistance(vox, i, whack);
		double curr = currentDistance(vox, i);
		double delta_dens = dDistanceTodDensity(AtomPtr(), curr, delta_dist);
		
		cuml += delta_dens;
	}
	
	return cuml;
}

double Gradiator::gradForWhack(WhackPtr whack)
{
	double cuml = 0;
	double cumlabs = 0;

	for (size_t i = 0; i < _voxels.size(); i++)
	{
		double add = deltaVoxel4Whack(whack, &_voxels[i]);
		cuml += add;
		cumlabs += fabs(add);
	}
	std::cout << "Cuml: " << cuml << std::endl;
	std::cout << "Cumlabs: " << cumlabs << std::endl;
	std::cout << "Net movement is " << cuml / cumlabs * 100
	<< "% of total movements" << std::endl;
	
	return cuml;
}

