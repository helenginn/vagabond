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

	double cutoff_d = 2.0;
	double cutoff_dsq = cutoff_d * cutoff_d;
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
	
	double cutoff = MAP_VALUE_CUTOFF;
	for (size_t j = 0; j < fft->nn; j++)
	{
		if (fft->data[j][1] < cutoff)
		{
			continue;
		}

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
			
			if (fabs(diff.x) > cutoff_d 
			    || fabs(diff.y) > cutoff_d 
			    || fabs(diff.z) > cutoff_d)
			{
				continue;
			}

			if (vec3_sqlength(diff) < cutoff_dsq)
			{
				ExplicitModelPtr m = a->getExplicitModel();
				std::vector<BondSample> *samples = m->getManyPositions();
				
				nearby.reserve(nearby.size() + samples->size());
				for (int k = 0; k < samples->size(); k++)
				{
					SingleAtom pair;
					pair.res = a->getResidueNum();
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
	
	double correl = correlationCoefficient();
	std::cout << "Current CC: " << correl << std::endl;
	
	for (int i = 0; i < _ws.size(); i++)
	{
		double nDelta = deltaCC(_ws[i], -1);
		double pDelta = deltaCC(_ws[i], 1);
		std::cout << i << ", " 
		<< std::setprecision(8) << nDelta << ", "
		<< std::setprecision(8) << pDelta << std::endl;
	}
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

double Gradiator::deltaVoxel4Whack(WhackPtr whack, Voxel *vox, int dir)
{
	std::vector<SingleAtom> *nearby = &vox->nearby_atoms;
	int whackRes = whack->getBond()->getAtom()->getResidueNum();
	
	double cuml = 0;
	int anchor = _polymer->getAnchor();
	bool left = whackRes < anchor;
	
	if (left) dir *= -1;
	
	for (size_t i = 0; i < nearby->size(); i++)
	{
		if (dir > 0)
		{
			if (whackRes >= nearby->at(i).res)
			{
				continue;
			}
		}
		else
		{
			if (whackRes < nearby->at(i).res)
			{
				continue;
			}
		}
		
		double delta_dist = deltaDistance(vox, i, whack);
		double curr = currentDistance(vox, i);
		double delta_dens = dDistanceTodDensity(AtomPtr(), curr, delta_dist);
		
		cuml += delta_dens;
	}
	
	return cuml;
}

double Gradiator::sum_x()
{
	double mean = 0;
	
	for (int i = 0; i < _voxels.size(); i++)
	{
		double calc = _voxels[i].calc;
		mean += calc;
	}
	
	mean /= (double)_voxels.size();
	
	return mean;
}

double Gradiator::sum_y()
{
	double mean = 0;
	
	for (int i = 0; i < _voxels.size(); i++)
	{
		double obs = _voxels[i].obs;
		mean += obs;
	}
	
	mean /= (double)_voxels.size();
	
	return mean;
}

double Gradiator::sxx()
{
	double sum = sum_x();
	double total = 0;
	
	for (int i = 0; i < _voxels.size(); i++)
	{
		double calc = _voxels[i].calc;
		double add = (calc - sum) * (calc - sum);
		total += add;
	}
	
	total /= (double)_voxels.size();
	return total;
}

double Gradiator::syy()
{
	double sum = sum_y();
	double total = 0;
	
	for (int i = 0; i < _voxels.size(); i++)
	{
		double obs = _voxels[i].obs;
		double add = (obs - sum) * (obs - sum);
		total += add;
	}
	
	total /= (double)_voxels.size();
	return total;
}

double Gradiator::sxy()
{
	double sy = sum_y();
	double sx = sum_x();

	double total = 0;
	
	for (int i = 0; i < _voxels.size(); i++)
	{
		double obs = _voxels[i].obs;
		double calc = _voxels[i].calc;
		double add = (obs - sy) * (calc - sx);
		total += add;
	}
	
	total /= (double)_voxels.size();
	return total;
}

double Gradiator::correlationCoefficient()
{
	double cc = sxy();
	double denom = sqrt(sxx() * syy());
	
	cc /= denom;
	
	return cc;
}

double Gradiator::deltaSumX4Whack(WhackPtr w, int dir)
{
	double sum = 0;
	
	for (int i = 0; i < _voxels.size(); i++)
	{
		double grad = deltaVoxel4Whack(w, &_voxels[i], dir);
		sum += grad;
	}
	
	sum /= (double)_voxels.size();
	return sum;
}

void Gradiator::deltaSs4Whack(WhackPtr w, double *dsxx, double *dsxy, int dir)
{
	double sx = sum_x();
	double sy = sum_y();
	double sum = deltaSumX4Whack(w, dir);

	*dsxx = 0;
	*dsxy = 0;

	for (int i = 0; i < _voxels.size(); i++)
	{
		double grad = deltaVoxel4Whack(w, &_voxels[i], dir);

		double add = 2 * (_voxels[i].calc - sx) * (grad - sum);
		*dsxx += add;

		add = (grad - sum) * (_voxels[i].obs - sy);
		*dsxy += add;
	}
	
	*dsxx /= (double)_voxels.size();
	*dsxy /= (double)_voxels.size();
}

double Gradiator::deltaCC(WhackPtr w, int dir)
{
	double deltaCC = 0;
	double xx = sxx();
	double denom_sq = xx * syy();
	
	double dsxy = 0; double dsxx = 0;
	deltaSs4Whack(w, &dsxx, &dsxy, dir);
	
	double left = dsxy * sqrt(xx);
	double right = sxy() * 0.5 / sqrt(xx) * dsxx;
	
	deltaCC = left - right;

	deltaCC *= sqrt(syy());
	deltaCC /= denom_sq;
	
	
	return deltaCC;
}

