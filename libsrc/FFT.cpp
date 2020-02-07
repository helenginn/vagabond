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

#include "FFT.h"
#include "Element.h"
#include "Atom.h"
#include "Timer.h"
#include "Anisotropicator.h"
#include "ExplicitModel.h"
#include <cstring>
#include <float.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>

#include "../libccp4/cmtzlib.h"
#include "../libccp4/csymlib.h"
#include "../libccp4/ccp4_spg.h"
#include "../libccp4/ccp4_general.h"

std::vector<FFTDim *> VagFFT::_dimensions;

VagFFT::VagFFT(VagFFT &fft)
{
	_nx = fft._nx;
	_ny = fft._ny;
	_nz = fft._nz;
	_nele = fft._nele;
	_nscratch = fft._nscratch;
	_activeScratch = fft._activeScratch;
	_nn = fft._nn;
	_status = fft._status;
	_stride = fft._stride;
	_total = fft._total;
	_origin = fft._origin;
	_toReal = fft._toReal;
	_toRecip = fft._toRecip;
	_recipBasis = fft._recipBasis;
	_realBasis = fft._realBasis;
	_setMatrices = fft._setMatrices;
	_myDims = fft._myDims;
	_elements = fft._elements;
	_elementMap = fft._elementMap;
	_spg = fft._spg;

	_data = (fftwf_complex *)fftwf_malloc(_total);
	int start = finalIndex(0);
	_lastData = &_data[start];

	if (!_data)
	{
		printf("ERROR: Malloc failed for VagFFT, nn = %i\n", _nn);
		exit(1);
	}

	memcpy(_data, fft._data, _total);
}

VagFFT::VagFFT(int nx, int ny, int nz, int nele, int scratches)
{
	_nx = nx;
	_ny = ny;
	_nz = nz;
	_myDims = NULL;
	_spg = CSym::ccp4spg_load_by_ccp4_num(1);

	if (_nx % 2 == 1) _nx -= 1;
	if (_ny % 2 == 1) _ny -= 1;
	if (_nz % 2 == 1) _nz -= 1;

	_nele = nele;
	_nscratch = scratches;
	_activeScratch = 0;
	_nn = _nx * _ny * _nz;
	_status = FFTEmpty;

	_stride = 2 * (_nele) + 1 + _nscratch;
	_total = _nn * _stride * sizeof(fftwf_complex);
	_data = (fftwf_complex *)fftwf_malloc(_total);

	int start = finalIndex(0);
	_lastData = &_data[start];

	if (!_data)
	{
		printf("ERROR: Malloc failed for VagFFT, nn = %i\n", _nn);
		exit(1);
	}

	memset(_data, 0, _total);

	_origin = empty_vec3();
	_toReal = make_mat3x3();
	_toRecip = make_mat3x3();
	_realBasis = make_mat3x3();
	_recipBasis = make_mat3x3();
	_setMatrices = false;
}

void VagFFT::addElement(ElementPtr ele)
{
	_elementMap[ele] = _elements.size();
	_elements.push_back(ele);
}

void VagFFT::printStatus()
{
	std::cout << "Current status: " << std::flush;
	
	switch (_status)
	{
		case FFTEmpty:
		std::cout << "FFTEmpty";
		break;
		case FFTRealSpace:
		std::cout << "FFTRealSpace";
		break;
		case FFTReciprocalSpace:
		std::cout << "FFTReciprocalSpace";
		break;
		case FFTSeparateAtoms:
		std::cout << "FFTSeperateAtoms";
		break;
		default:
		std::cout << "Unknown";
		break;
	}
	
	std::cout << std::endl;
}

bool VagFFT::sanityCheck()
{
	if (_nx <= 0 || _ny <= 0 || _nz <= 0)
	{
		std::cout << "Illogical VagFFT dimensions: " << _nx 
		<< " " << _ny << " " << _nz << std::endl;
		return false;
	}
	
	if (_nele != _elements.size())
	{
		std::cout << "Element number and assigned element count do"\
		" not match, " << _nele << " vs. " << _elements.size() << std::endl;
		return false;
	}
	
	if (!_setMatrices)
	{
		std::cout << "Need to set up real / reciprocal change of bases" 
		<< std::endl;
		return false;
	}
	
	return true;
}

void VagFFT::wipe()
{
	memset(_data, 0, _total);
	_status = FFTEmpty;
}

void VagFFT::makePlans()
{
	if (!sanityCheck())
	{
		std::cout << "Sanity check failed" << std::endl;
		return;
	}

	/* do we already have this planned dimension? */
	for (int i = 0; i < _dimensions.size(); i++)
	{
		if (_dimensions[i]->nx == _nx
		    && _dimensions[i]->ny == _ny
		    && _dimensions[i]->nz == _nz
		    && _dimensions[i]->nele == _nele)
		{
			_myDims = _dimensions[i];
			return;
		}
	}
	
	Timer timer;
	if (_nele > 0)
	{
		std::cout << "Planning FFT... " << std::flush;
	}


	FFTDim *dims = (FFTDim *)malloc(sizeof(FFTDim));
	dims->nx = _nx; dims->ny = _ny; dims->nz = _nz; dims->nele = _nele;
	_dimensions.push_back(dims);
	_myDims = _dimensions.back();
	
	int ns[3] = {_nz, _ny, _nx};
	unsigned fftw_flags = FFTW_MEASURE;

	fftwf_plan plan = fftwf_plan_many_dft(3, ns, _nele, _data, NULL, 
	                                      _stride, 2, _data, 
	                                      NULL, _stride, 2, 
	                                      1, fftw_flags); 
	_myDims->atom_real_to_recip = plan;

	plan = fftwf_plan_many_dft(3, ns, _nele, _data, NULL, _stride, 
	                                      2, _data, 
	                                      NULL, _stride, 2, 
	                                      -1, fftw_flags); 
	_myDims->atom_recip_to_real = plan;

	plan = fftwf_plan_many_dft(3, ns, 1, _lastData, NULL, /*null=embed*/
	                           _stride, 2, _lastData, NULL, 
	                           _stride, 2, 1, fftw_flags); 

	_myDims->real_to_recip = plan;

	plan = fftwf_plan_many_dft(3, ns, 1, _lastData, NULL, 
	                           _stride, 2, _lastData, NULL,
	                           _stride, 2, -1, fftw_flags); 

	_myDims->recip_to_real = plan;
	
	if (_nele > 0)
	{
		std::cout << " done. ";
		timer.quickReport();
		std::cout << std::endl;
	}
}

void VagFFT::multiplyDotty(float val)
{
	for (int i = 0; i < _nn; i++)
	{
		for (int j = 0; j < _nele; j++)
		{
			long dotty_index = dottyIndex(i, j);
			_data[dotty_index][0] *= val;
			_data[dotty_index][1] *= val;
		}
	}
}

void VagFFT::multiplyFinal(float val)
{
	for (int i = 0; i < _nn; i++)
	{
		long final_index = finalIndex(i);
		_data[final_index][0] *= val;
		_data[final_index][1] *= val;
	}
}


void VagFFT::setupElements(bool wipe)
{
	if (_elements.size() == 0 && _nele == 0)
	{
		return;
	}

	for (int z = -_nz / 2; z <= _nz / 2; z++)
	{
		for (int y = -_ny / 2; y <= _ny / 2; y++)
		{
			for (int x = -_nx / 2; x <= _nx / 2; x++)
			{
				long ele = element(x, y, z);
				
				vec3 real = make_vec3(x, y, z);
				mat3x3_mult_vec(_toRecip, &real);
				mat3x3_mult_vec(_realBasis, &real);

				for (int j = 0; j < _nele; j++)
				{
					long ele_index = eleIndex(ele, j);
					
					ElementPtr elem = _elements[j];
					double val = Element::getVoxelValue(&*elem, real.x, 
					                                    real.y, real.z);

					_data[ele_index][0] = val;
					_data[ele_index][1] = 0;
					
					if (wipe)
					{
						long dotty_index = ele_index - 1;
						_data[dotty_index][0] = 0; 
						_data[dotty_index][1] = 0;
					}
				}
			}
		}
	}
}

void VagFFT::prepareAtomSpace()
{
	if (_status == FFTEmpty)
	{
		setupElements();
	}
	else
	{
		for (int i = 0; i < _nn; i++)
		{
			// elements should be hoovered up straight after last

			// transform.

			long final_index = finalIndex(i);
			_data[final_index][0] = 0; 
			_data[final_index][1] = 0;
		}
	}

	_status = FFTSeparateAtoms;
}

void VagFFT::separateAtomTransform()
{
	fftwf_execute_dft(_myDims->atom_real_to_recip, _data, _data);
	
	/* go through and multiply transform of atoms by element factor */
	for (int i = 0; i < _nn; i++)
	{
		for (int j = 0; j < _nele; j++)
		{
			long dotty_index = dottyIndex(i, j);
			long ele_index = dotty_index + 1;

			_data[dotty_index][0] *= _data[ele_index][0];
			_data[dotty_index][1] *= _data[ele_index][0];
		}
	}

	/* add up final indices as sum of individual atom types */
	for (int i = 0; i < _nn; i++)
	{
		long final_index = finalIndex(i);
		_data[final_index][0] = 0;
		_data[final_index][1] = 0;

		for (int j = 0; j < _nele; j++)
		{
			long dotty_index = dottyIndex(i, j);

			_data[final_index][0] += _data[dotty_index][0];
			_data[final_index][1] += _data[dotty_index][1];
			
			_data[dotty_index][0] = 0;
			_data[dotty_index][1] = 0;
		}
	}

	_status = FFTReciprocalSpace;
}

void VagFFT::fft(FFTTransform transform)
{
	if (_myDims == NULL)
	{
		std::cout << "Attempting to transform FFT but not set up any"
		" plans for doing so." << std::endl;
		exit(1);
	}

	if (transform == FFTAtomsToReciprocal || transform == FFTAtomsToReal)
	{
		if (!(_status == FFTSeparateAtoms || _status == FFTEmpty))
		{
			std::cout << "Trying to transform individual atoms when status is "
			<< "not FFTSeparateAtoms" << std::endl;
			printStatus();
			exit(1);
		}

		separateAtomTransform();

		if (transform == FFTAtomsToReal)
		{
			fftwf_execute_dft(_myDims->recip_to_real, _lastData, _lastData);
			multiplyFinal(1 / (double)_nn);
			_status = FFTRealSpace;
		}
	}

	if (transform == FFTRealToReciprocal)
	{
		if (!(_status == FFTRealSpace || _status == FFTEmpty))
		{
			std::cout << "Asking for transform Real to Reciprocal, but "\
			"FFT status is not Real Space" << std::endl;
			printStatus();
			exit(1);
		}

		fftwf_execute_dft(_myDims->real_to_recip, _lastData, _lastData);
		_status = FFTReciprocalSpace;
	}

	if (transform == FFTReciprocalToReal)
	{
		if (!(_status == FFTReciprocalSpace || _status == FFTEmpty))
		{
			std::cout << "Asking for transform Reciprocal to Real, but "\
			"FFT status is not Reciprocal Space" << std::endl;
			printStatus();
			exit(1);
		}

		fftwf_execute_dft(_myDims->recip_to_real, _lastData, _lastData);
		multiplyFinal(1 / (double)_nn);
		_status = FFTRealSpace;
	}

}

void VagFFT::setUnitCell(std::vector<double> dims)
{
	_toReal = mat3x3_from_unit_cell(dims[0], dims[1], dims[2], dims[3],
	                                dims[4], dims[5]);
	_toRecip = mat3x3_inverse(_toReal);

	_realBasis = _toReal;
	mat3x3_scale(&_realBasis, 1 / (double)_nx, 1 / (double)_ny, 
	             1 / (double)_nz);
	_recipBasis = mat3x3_inverse(_realBasis);
	_setMatrices = true;
}

void VagFFT::addImplicitAtom(AtomPtr atom)
{
	mat3x3 tensor = atom->getModel()->getRealSpaceTensor();

	Anisotropicator aniso;
	aniso.setTensor(tensor);
	mat3x3 ellipsoid = aniso.basis();
	
	/* make these variances into standard deviations */
	for (int i = 0; i < 3; i++)
	{
		vec3 axis = mat3x3_axis(ellipsoid, i);
		double length = vec3_length(axis);
		
		for (int j = 0; j < 3; j++)
		{
			ellipsoid.vals[j * 3 + i] *= 1 / sqrt(length);
		}
	}

	vec3 centre = atom->getAbsolutePosition();
	
	/* TODO find limiting voxels */
	const int scale = 5.0; /* no. standard devs we care about */

	vec3 maxVals = make_vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	
	for (int i = 0; i < 3; i++)
	{
		vec3 axis = mat3x3_axis(ellipsoid, i);
		double length = vec3_length(axis);
		*(&maxVals.x + i) = length;
	}
	
	vec3_mult(&maxVals, scale);
	mat3x3_mult_vec(_recipBasis, &maxVals);
	vec3 four = make_vec3(4, 4, 4);
	vec3_max_each(&maxVals, four);
	
	vec3_subtract_from_vec3(&centre, _origin);
	mat3x3_mult_vec(_recipBasis, &centre);
	ElementPtr ele = atom->getElement();

	double occ = atom->getModel()->getEffectiveOccupancy();
	double total = populateImplicit(ele, centre, maxVals, 
	                                ellipsoid, occ, true);
}

double VagFFT::populateImplicit(ElementPtr ele, vec3 centre, vec3 maxVals,
                                mat3x3 tensor, double scale, bool add)
{
	double total = 0;
	double count = 0;
	
	mat3x3 inv = mat3x3_inverse(tensor);
	vec3 xl = mat3x3_axis(tensor, 0);
	vec3 yl = mat3x3_axis(tensor, 1);
	vec3 zl = mat3x3_axis(tensor, 2);
	double std_x = vec3_length(xl);
	double std_y = vec3_length(yl);
	double std_z = vec3_length(zl);
	double vol = mat3x3_volume(_realBasis);
	int column = whichColumn(ele);
	
	for (int z = centre.z - maxVals.z; z <= centre.z + maxVals.z; z++)
	{
		for (int y = centre.y - maxVals.y; y <= centre.y + maxVals.y; y++)
		{
			for (int x = centre.x - maxVals.x; x <= centre.x + maxVals.x; x++)
			{
				vec3 pos = make_vec3(x, y, z);
				/* in voxel system */
				vec3 offset = vec3_subtract_vec3(pos, centre);
				/* and then to Angstroms */
				mat3x3_mult_vec(_realBasis, &offset);
				/* to aniso-U-sensitive coordinates */
				mat3x3_mult_vec(inv, &offset);

				offset.x = fabs(offset.x);
				offset.y = fabs(offset.y);
				offset.z = fabs(offset.z);
				offset.x *= offset.x / 2;
				offset.y *= offset.y / 2;
				offset.z *= offset.z / 2;
				
				double factor = 1 / sqrt(2 * M_PI);
				double dens = (factor / std_x) * exp(- offset.x);
				dens *= (factor / std_y) * exp(- offset.y);
				dens *= (factor / std_z) * exp(- offset.z);

				if (add)
				{
					addInterpolatedToReal(column, x, y, z, dens * scale);
				}
				
				total += dens * vol;
			}
		}
	}

	return total;
}

int VagFFT::whichColumn(ElementPtr ele)
{
	if (_elementMap.count(ele))
	{
		return _elementMap[ele] * 2;
	}

	return -1;
}

void VagFFT::addInterpolatedToReal(int column, double sx, double sy, 
                                   double sz, double val)
{
	int lx = (int)floor(sx);
	int ly = (int)floor(sy);
	int lz = (int)floor(sz);

	double xProps[2];
	double yProps[2];
	double zProps[2];

	xProps[1] = sx - (double)lx;
	yProps[1] = sy - (double)ly;
	zProps[1] = sz - (double)lz;

	xProps[0] = 1 - xProps[1];
	yProps[0] = 1 - yProps[1];
	zProps[0] = 1 - zProps[1];
	
	for (int r = 0; r < 2; r++)
	{
		for (int q = 0; q < 2; q++)
		{
			for (int p = 0; p < 2; p++)
			{
				int sx1 = lx + p;
				int sy1 = ly + q;
				int sz1 = lz + r;

				long index = element(sx1, sy1, sz1);
				double prop = xProps[p] * yProps[q] * zProps[r];
				_data[index * _stride + column][0] += prop * val;
			}	
		}
	}
}

double VagFFT::getPhase(int x, int y, int z)
{
	long pre_index = element(x, y, z);
	long index = finalIndex(pre_index);

	double degrees = atan2(_data[index][1], _data[index][0]) * 180 / M_PI;

	while (degrees >= 360) degrees -= 360;

	while (degrees < 0) degrees += 360;

	return degrees;

}

double VagFFT::getAmplitude(ElementPtr ele, int x, int y, int z)
{
	int col = whichColumn(ele);
	long index = element(x, y, z) * _stride + col + 1;

	double val =  (_data[index][0] * _data[index][0] + 
	               _data[index][1] * _data[index][1]);
	
	return sqrt(val);
}

void VagFFT::printSlice(double zVal, double scale)
{
	sanityCheck();

	if (scale < 0)
	{
		scale = averageAll() * 5;
	}
	
	zVal *= _nz;
	
	for (int j = 0; j < _ny; j++)
	{
		std::cout << "| ";
		for (int i = 0; i < _nx; i++)
		{
			std::string symbol = " ";
			double value = getAmplitude(_elements[0], i, j, zVal);

			if (value > 0.01 * scale) symbol = ".";
			if (value > 0.02 * scale) symbol = ":";
			if (value > 0.04 * scale) symbol = "\"";
			if (value > 0.08 * scale) symbol = "*";
			if (value > 0.16 * scale) symbol = "x";
			if (value > 0.32 * scale) symbol = "H";
			if (value > 0.64 * scale) symbol = "#";
			if (value > 1.00 * scale) symbol = "@";

			std::cout << symbol;
		}

		std::cout << " |" << std::endl;
	}
	std::cout << std::endl;

}

void VagFFT::addExplicitAtom(AtomPtr atom)
{
	std::vector<BondSample> positions;
	positions = atom->getExplicitModel()->getFinalPositions();
	ElementPtr ele = atom->getElement();
	int column = whichColumn(ele);
	double vol = mat3x3_volume(_realBasis);
	double low = atom->getExplicitModel()->getLowestZ();

	for (int i = 0; i < positions.size(); i++)
	{
		vec3 pos = positions[i].start;

		/* force cache hit for the rest, hopefully */
		if (i == 0)
		{
			vec3 first = make_vec3(pos.x - _origin.x, 
			                       pos.y - _origin.y, 
			                       low - _origin.z);

			mat3x3_mult_vec(_recipBasis, &first);
			addInterpolatedToReal(column, first.x, first.y, first.z, 0.);
		}

		vec3_subtract_from_vec3(&pos, _origin);
		mat3x3_mult_vec(_recipBasis, &pos);

		double occ = positions[i].occupancy;
		double dens = occ / vol;

		addInterpolatedToReal(column, pos.x, pos.y, pos.z, dens);
	}
}

void VagFFT::addAtom(AtomPtr atom)
{
	if (atom->getModel()->hasExplicitPositions())
	{
		addExplicitAtom(atom);
	}
	else
	{
		addImplicitAtom(atom);
	}
}

/* 11-point interpolation - attempted transcription from Dave's function
 * from GAP */
double VagFFT::cubic_interpolate(vec3 vox000, size_t im)
{
	/* vox000 has integer and real components */
	
	/* Pick out just the real components - this is faster
	 * than modf */
	vec3 uvw = make_vec3(vox000.x - (double)((int)vox000.x),
	                        vox000.y - (double)((int)vox000.y),
	                        vox000.z - (double)((int)vox000.z));

	/* Extra refers to the additional index to be fished
	 * for 11-point interpolation. We already get 0 and 1. */
	int extra[3] = {-1, -1, -1};
	int central[3] = {0, 0, 0};
	int next[3] = {1, 1, 1};
	
	/* If uvw components are greater than 0.5, then flip them 
	 * make the extra index one ahead and reverse the order */
	for (int i = 0; i < 3; i++)
	{
		if (*(&uvw.x + i) > 0.5)
		{
			extra[i] = 2;
			central[i] = 1;
			next[i] = 0;
			*(&uvw.x + i) = 1 - *(&uvw.x + i);
		}
	}

	int vox000x = vox000.x + central[0];
	int vox000y = vox000.y + central[1];
	int vox000z = vox000.z + central[2];
	int vox000xm = vox000.x + next[0];
	int vox000ym = vox000.y + next[1];
	int vox000zm = vox000.z + next[2];
	int vox000xn = vox000.x + extra[0];
	int vox000yn = vox000.y + extra[1];
	int vox000zn = vox000.z + extra[2];

	collapse(&vox000x, &vox000y, &vox000z);
	collapse(&vox000xm, &vox000ym, &vox000zm);
	collapse(&vox000xn, &vox000yn, &vox000zn);

	vox000y  *= _nx;
	vox000ym *= _nx;
	vox000yn *= _nx;
	vox000z  *= _nx * _ny;
	vox000zm *= _nx * _ny;
	vox000zn *= _nx * _ny;

	long idx000 = vox000x + vox000y + vox000z;
	long idx100 = vox000xm + vox000y + vox000z;
	long idx010 = vox000x + vox000ym + vox000z;
	long idx110 = vox000xm + vox000ym + vox000z;
	long idx001 = vox000x + vox000y + vox000zm;
	long idx101 = vox000xm + vox000y + vox000zm;
	long idx011 = vox000x + vox000ym + vox000zm;
	long idx111 = vox000xm + vox000ym + vox000zm;
	
	long idxn00 = vox000xn + vox000y + vox000z;
	long idx0n0 = vox000x + vox000yn + vox000z;
	long idx00n = vox000x + vox000y + vox000zn;
	
	double u = uvw.x;
	double v = uvw.y;
	double w = uvw.z;
	
	double p000 = _data[finalIndex(idx000)][im];
	double p001 = _data[finalIndex(idx001)][im];
	double p010 = _data[finalIndex(idx010)][im];
	double p011 = _data[finalIndex(idx011)][im];
	double p100 = _data[finalIndex(idx100)][im];
	double p101 = _data[finalIndex(idx101)][im];
	double p110 = _data[finalIndex(idx110)][im];
	double p111 = _data[finalIndex(idx111)][im];
	
	double a = p100 - p000;
	double b = p010 - p000;
	double c = p110 - p010;
	double d = p101 - p001;
	
	double pn00 = _data[finalIndex(idxn00)][im];
	double p0n0 = _data[finalIndex(idx0n0)][im];
	double p00n = _data[finalIndex(idx00n)][im];

	double p8value = p000+u*(a+w*(-a+d)+v*((c-a)+w*( a-c-d-p011+p111)))
	+ v*(b+w*(-p001+p011-b))+w*(-p000+p001);
	
	double mod = (p000 - 0.5 * p100 - 0.5 * pn00) * (u - u * u);
	mod += (p000 - 0.5 * p010 - 0.5 * p0n0) * (v - v * v);
	mod += (p000 - 0.5 * p001 - 0.5 * p00n) * (w - w * w);
	
	double p11value = p8value + 0.4 * mod;

	return p11value;
}

/*  For multiplying point-wise
 *  No assumption that interpolation is not needed.
 */
double VagFFT::operation(VagFFTPtr fftCrystal, VagFFTPtr fftAtom,
                      MapScoreType mapScoreType, std::vector<CoordVal> *vals,
                      bool sameScale, bool fcOnly)
{
	/* I rarely comment something so heavily but I will get confused if
	 * I don't, this time, as I can't soak the protocol into the variable
	 * names. Bear in mind the three coordinate systems:
	 * (a) Angstroms
	 * (b) Crystal voxels
	 * (c) Atom voxels */

	/* find offset from origins of two maps */
	vec3 add = vec3_subtract_vec3(fftAtom->_origin, fftCrystal->_origin);

	/* find the centre of the atom FFT grid */
	vec3 half_box = make_vec3(0.5, 0.5, 0.5);
	mat3x3_mult_vec(fftAtom->_toReal, &half_box);
	vec3_add_to_vec3(&add, half_box);

	/* to reciprocal space coordinates */
	mat3x3_mult_vec(fftCrystal->_toRecip, &add);

	/* Bring the fractional coordinate of the atom into range 0 < frac <= 1 */
	FFT::collapseFrac(&add.x, &add.y, &add.z);

	/* Multiply by the relative dimensions of the crystal */
	double multX = add.x * fftCrystal->_nx;
	double multY = add.y * fftCrystal->_ny;
	double multZ = add.z * fftCrystal->_nz;

	/* Get the remainder after subtracting a number of whole voxels */
	vec3 crystFracOffset; // in crystal coordinates, converted later.
	crystFracOffset.x = fmod(multX, 1.);
	crystFracOffset.y = fmod(multY, 1.);
	crystFracOffset.z = fmod(multZ, 1.);

	/* Store the non-remainder whole voxel values for way later. */
	vec3 atomWholeCrystCoords = make_vec3((int)multX, (int)multY, (int)multZ);

	/* Prepare a matrix to convert crystal voxels into atomic voxels */
	mat3x3 crystal2AtomVox = mat3x3_mult_mat3x3(fftAtom->getRecipBasis(),
	                                            fftCrystal->getRealBasis());

	/* Prepare a matrix to convert atomic voxels into crystal voxels */
	mat3x3 atomVox2Crystal = mat3x3_mult_mat3x3(fftCrystal->getRecipBasis(),
	                                            fftAtom->getRealBasis());

	/* Apply this offset and reverse it. 
	 * This small offset must be added to all future atomic coordinates prior
	 * to interpolation. This will be subtracted from an atom voxel to find
	 * the relevant crystal voxel, and should therefore be in atom voxels */
	mat3x3_mult_vec(crystal2AtomVox, &crystFracOffset);
	vec3_mult(&crystFracOffset, -1);

	/* There will be an additional shift having moved the atom by
	 * half the dimension length which needs to be taken into account, 
	 * unfortunately. */
	vec3 atomShift = make_vec3((double)(-fftAtom->_nx) * 0.5,
	                           (double)(-fftAtom->_ny) * 0.5,
	                           (double)(-fftAtom->_nz) * 0.5);
	vec3 shift = mat3x3_mult_vec(atomVox2Crystal, atomShift);

	/* In crystal voxels at the moment - don't worry about fractional
	 * shifts. */
	vec3 shiftRemainder = make_vec3(fmod(shift.x, 1),
	                                fmod(shift.y, 1),
	                                fmod(shift.z, 1));

	vec3 wholeShiftOnly = make_vec3(shift.x - shiftRemainder.x - 1,
	                                shift.y - shiftRemainder.y - 1,
	                                shift.z - shiftRemainder.z - 1);

	shiftRemainder.x = 1 + shiftRemainder.x;
	shiftRemainder.y = 1 + shiftRemainder.y;
	shiftRemainder.z = 1 + shiftRemainder.z;

	/* The crystal voxels must be converted to atomic voxels to determine
	 * final offset for atom sampling. */
	mat3x3_mult_vec(crystal2AtomVox, &shiftRemainder);
	vec3_mult(&shiftRemainder, -1);

	/* Corner of atom as whole value offset in the crystal coordinates. */
	vec3 cornerCrystal = vec3_add_vec3(wholeShiftOnly, atomWholeCrystCoords);
//	vec3 cornerCrystal = atomWholeCrystCoords;

	/* Fractional offset in atomic coordinates, for each atom as a
	 * 	fraction of the crystal voxel. */
	crystFracOffset = vec3_add_vec3(crystFracOffset, shiftRemainder);
	vec3 crystOffset = mat3x3_mult_vec(atomVox2Crystal, crystFracOffset);
	
	/* We loop around these crystal voxel limits now (ss -> ms -> fs).
	 * We also break the loop if it exceeds the limits of our atom voxels
	 * during the loop itself. */

	/* Determine bounding box - 9th Dec 2017 */
	vec3 minAtom = make_vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	vec3 maxAtom = make_vec3(0, 0, 0);

	for (int k = 0; k <= fftAtom->_nz; k += fftAtom->_nz)
	{
		for (int j = 0; j <= fftAtom->_ny; j += fftAtom->_ny)
		{
			for (int i = 0; i <= fftAtom->_nx; i += fftAtom->_nx)
			{
				vec3 atomCorner = make_vec3(i, j, k);
				vec3 toCrystal = mat3x3_mult_vec(atomVox2Crystal, atomCorner);
				toCrystal.x = (int)toCrystal.x;
				toCrystal.y = (int)toCrystal.y;
				toCrystal.z = (int)toCrystal.z;

				vec3_min_each(&minAtom, toCrystal);
				vec3_max_each(&maxAtom, toCrystal);
			}
		}
	}
	
	double step = 1;
	long count = 0;

	/* min/maxAtoms are in crystal coordinates.*/
	for (double k = minAtom.z; k < maxAtom.z; k += step)
	{
		for (double j = minAtom.y; j < maxAtom.y; j += step)
		{
			for (double i = minAtom.x; i < maxAtom.x; i += step)
			{
				/* Position currently in crystal coords - change to atom. */
				vec3 crystalPos = make_vec3(i, j, k);
				vec3 atomPos = crystalPos;

				if (!sameScale)
				{
					mat3x3_mult_vec(crystal2AtomVox, &atomPos);

					if (atomPos.x < 0 || atomPos.y < 0 || atomPos.z < 0)
					{
						continue;
					}

					if (atomPos.x > fftAtom->_nx || atomPos.y > fftAtom->_ny
					    || atomPos.z > fftAtom->_nz)
					{
						continue;
					}
				}

				/* Now we must find the relative crystal voxel to write this
				 * density value to, given that the atom was wrapped around
				 * the origin (center). This should work regardless of odd/
				 * even dimension lengths. */

				/* We add the tiny offset which resulted from the atom
				 * falling between two voxels, in atomic voxels */
				vec3_add_to_vec3(&atomPos, crystFracOffset);
				
				/* Find the interpolated value which atomPos falls on */
				double atomReal = 0; double atomImag = 0;

				while (atomPos.x < 0) atomPos.x += fftAtom->_nx;
				while (atomPos.y < 0) atomPos.y += fftAtom->_ny;
				while (atomPos.z < 0) atomPos.z += fftAtom->_nz;

				while (atomPos.x >= fftAtom->_nx) atomPos.x -= fftAtom->_nx;
				while (atomPos.y >= fftAtom->_ny) atomPos.y -= fftAtom->_ny;
				while (atomPos.z >= fftAtom->_nz) atomPos.z -= fftAtom->_nz;

				atomReal = fftAtom->cubic_interpolate(atomPos, 0);

				if (vals != NULL)
				{
					atomImag = fftAtom->cubic_interpolate(atomPos, 1);
				}

				/* We add the crystal offset so we don't end up with thousands
				 * of atoms at the very centre of our map */
				vec3 cVox = empty_vec3();
				long cIndex = 0;
				
				if (!fcOnly)
				{
					cVox = vec3_add_vec3(crystalPos, cornerCrystal);

					if (mapScoreType == MapScoreAddNoWrap)
					{
						if (cVox.x < -fftCrystal->_nx / 2 || 
						    cVox.y < -fftCrystal->_ny / 2 ||
						    cVox.z < -fftCrystal->_nz / 2)
						{
							continue;
						}					

						if (cVox.x > fftCrystal->_nx / 2 ||
						    cVox.y > fftCrystal->_ny / 2 ||
						    cVox.z > fftCrystal->_nz / 2)
						{
							continue;
						}
					}

					while (cVox.x < 0) cVox.x += fftCrystal->_nx;
					while (cVox.y < 0) cVox.y += fftCrystal->_ny;
					while (cVox.z < 0) cVox.z += fftCrystal->_nz;

					while (cVox.x>=fftCrystal->_nx) cVox.x -= fftCrystal->_nx;
					while (cVox.y>=fftCrystal->_ny) cVox.y -= fftCrystal->_ny;
					while (cVox.z>=fftCrystal->_nz) cVox.z -= fftCrystal->_nz;

					/* Get the index of this final crystal voxel. */
					cIndex = fftCrystal->element(cVox.x + 0.5,
					                             cVox.y + 0.5,
					                             cVox.z + 0.5);
				}

				if (mapScoreType == MapScoreTypeCorrel)
				{
					if (!fcOnly)
					{
						// We do NOT need to interpolate //
						double realCryst = fftCrystal->getReal(cIndex);

						if (vals)
						{
							CoordVal val;
							val.fo = realCryst;
							val.fc = atomReal;
							val.weight = atomImag;
							#ifdef COORDVAL_FULL
							vec3 frac = fftCrystal->fracFromElement(cIndex);
							val.pos = frac;
							#endif

							vals->push_back(val);
						}
					}
					else
					{
						vals->at(count).fc = atomReal;
						vals->at(count).weight = atomImag;
						count++;
					}
				}
				else if (mapScoreType == MapScoreTypeNone)
				{
					/* Add the density to the real value of the crystal
					 * voxel. */

					fftCrystal->addToReal(cIndex, atomReal);
				}
			}
		}
	}

	return 0;
}

void VagFFT::setScale(double cubeDim)
{
	mat3x3 real = make_mat3x3();
	mat3x3_scale(&real, cubeDim, cubeDim, cubeDim);
	_realBasis = real;
	_recipBasis = mat3x3_inverse(real);
	
	_toReal = _realBasis;
	mat3x3_scale(&_toReal, _nx, _ny, _nz);

	_toRecip = mat3x3_inverse(_toReal);
	_setMatrices = true;
}

void VagFFT::addSimple(VagFFTPtr v2)
{
	for (int i = 0; i < _nn; i++)
	{
		long final_index = finalIndex(i);
		_data[final_index][0] += v2->_data[final_index][0];
	}
}

void VagFFT::addSimple(FFTPtr v2)
{
	for (int i = 0; i < _nn; i++)
	{
		long final_index = finalIndex(i);
		_data[final_index][0] += v2->data[i][0];
	}
}

void VagFFT::copyFrom(FFTPtr fft)
{
	wipe();
	for (int i = 0; i < _nn; i++)
	{
		long final_index = finalIndex(i);
		_data[final_index][0] += fft->data[i][0];
		_data[final_index][1] += fft->data[i][1];
	}
}

void VagFFT::copyRealToImaginary()
{
	for (int i = 0; i < _nn; i++)
	{
		long final_index = finalIndex(i);
		_data[final_index][1] = _data[final_index][0];
	}
}

double VagFFT::sumReal(int scratch)
{
	double sum = 0;

	for (int i = 0; i < _nn; i++)
	{
		long index = 0;
		if (scratch < 0 || scratch >= _nscratch)
		{
			index = finalIndex(i);
		}
		else
		{
			index = scratchIndex(i, scratch);
		}

		sum += fabs(_data[index][0]);
	}

	return sum;
}

double VagFFT::getAmplitude(long i)
{
	long index = finalIndex(i);
	double val = _data[index][0] * _data[index][0] + 
	_data[index][1] * _data[index][1];
	return sqrt(val);
}

vec3 VagFFT::fracFromElement(long int element)
{
	long x = element % _nx;
	element -= x;
	element /= _nx;

	long y = element % _ny;
	element -= y;
	element /= _ny;

	long z = element;

	double xfrac = (double)x / (double)_nx;
	double yfrac = (double)y / (double)_ny;
	double zfrac = (double)z / (double)_nz;

	return make_vec3(xfrac, yfrac, zfrac);
}

void VagFFT::writeToFile(std::string filename, double maxResolution,
                         FFTPtr data, VagFFTPtr diff, VagFFTPtr calc)
{
	if (_status != FFTReciprocalSpace)
	{
		std::cout << "Warning: writing MTZ file for grid in "
		"real space: " << filename << std::endl;
	}
	
	double nLimit[3];
	nLimit[0] = _nx;
	nLimit[1] = _ny;
	nLimit[2] = _nz;
	
	for (int i = 0; i < 3; i++)
	{
		nLimit[i] = nLimit[i] - ((int)nLimit[i] % 2); // make even
		nLimit[i] /= 2;
	}
	
	if (data)
	{
		// lower limit if data has smaller spacing
		nLimit[0] = (nLimit[0] > data->nx) ? data->nx : nLimit[0];
		nLimit[1] = (nLimit[1] > data->ny) ? data->ny : nLimit[1];
		nLimit[2] = (nLimit[2] > data->nz) ? data->nz : nLimit[2];
	}
	
	double dStar = 1 / maxResolution;

	if (maxResolution <= 0) dStar = FLT_MAX;

	/* For writing MTZ files */

	int columns = 12;

	float cell[6], wavelength;
	float *fdata = new float[columns];

	/* variables for symmetry */
	float rsm[192][4][4];
	char ltypex[2];

	/* variables for MTZ data structure */
	CMtz::MTZ *mtzout;
	CMtz::MTZXTAL *xtal;
	CMtz::MTZSET *set;
	CMtz::MTZCOL *colout[columns + 1];

	double unitCell[6];
	unit_cell_from_mat3x3(_toReal, unitCell);
	cell[0] = unitCell[0];
	cell[1] = unitCell[1];
	cell[2] = unitCell[2];
	cell[3] = unitCell[3];
	cell[4] = unitCell[4];
	cell[5] = unitCell[5];
	wavelength = 1.00; // fixme

	std::string outputFileOnly = filename;
	std::string outputFile = FileReader::addOutputDirectory(outputFileOnly);

	mtzout = CMtz::MtzMalloc(0, 0);
	ccp4_lwtitl(mtzout, "Written from Helen's XFEL tasks ", 0);
	mtzout->refs_in_memory = 0;
	mtzout->fileout = CMtz::MtzOpenForWrite(outputFile.c_str());

	if (!mtzout->fileout)
	{
		std::cout << "Could not open " << outputFile.c_str() << std::endl;
		std::cerr << "Error: " << strerror(errno) << std::endl;
		return;
	}

	// then add symm headers...
	for (int i = 0; i < _spg->nsymop; ++i)
	CCP4::rotandtrn_to_mat4(rsm[i], _spg->symop[i]);
	strncpy(ltypex, _spg->symbol_old, 1);
	ccp4_lwsymm(mtzout, _spg->nsymop, _spg->nsymop_prim, rsm, ltypex,
	            _spg->spg_ccp4_num, _spg->symbol_old, _spg->point_group);

	// then add xtals, datasets, cols
	xtal = MtzAddXtal(mtzout, "vagabond_crystal", "vagabond_project", cell);
	set = MtzAddDataset(mtzout, xtal, "Dataset", wavelength);
	colout[0] = MtzAddColumn(mtzout, set, "H", "H");
	colout[1] = MtzAddColumn(mtzout, set, "K", "H");
	colout[2] = MtzAddColumn(mtzout, set, "L", "H");
	colout[3] = MtzAddColumn(mtzout, set, "FREE", "I");
	colout[4] = MtzAddColumn(mtzout, set, "FP", "F");
	colout[5] = MtzAddColumn(mtzout, set, "SIGFP", "Q");
	colout[6] = MtzAddColumn(mtzout, set, "FC", "F");
	colout[7] = MtzAddColumn(mtzout, set, "FWT", "F");
	colout[8] = MtzAddColumn(mtzout, set, "PHIC", "P");
	colout[9] = MtzAddColumn(mtzout, set, "PHWT", "P");
	colout[10] = MtzAddColumn(mtzout, set, "DELFWT", "F");
	colout[11] = MtzAddColumn(mtzout, set, "PHDELWT", "P");

	int num = 0;

	/* symmetry issues */
	for (int k = -nLimit[2]; k < nLimit[2]; k++)
	{
		for (int j = -nLimit[1]; j < nLimit[1]; j++)
		{
			for (int i = -nLimit[0]; i < nLimit[0]; i++)
			{
				bool asu = CSym::ccp4spg_is_in_asu(_spg, i, j, k);
				bool f000 = (i == 0 && j == 0 && k == 0);
				int abs = CSym::ccp4spg_is_sysabs(_spg, i, j, k);
				
				if (!asu || f000 || abs)
				{
					continue;
				}

				vec3 pos = make_vec3(i, j, k);
				mat3x3_mult_vec(_toRecip, &pos);

				if (vec3_length(pos) > dStar)
				{
					continue;
				}

				/* weighted amplitude and phase */
				double fwt = getAmplitude(i, j, k);
				double phwt = getPhase(i, j, k);

				/* get calculated amplitude from calc map, ideally
				 * available */
				double calcAmp = getAmplitude(i, j, k);
				double phic = phwt;
				
				if (calc)
				{
					calcAmp = calc->getAmplitude(i, j, k);
					phic = calc->getPhase(i, j, k);
				}

				int free = 1;

				/* Sort out observed values from observed diffraction,
				 * ideally available */
				double foInt = fwt * fwt;
				double foAmp = fwt;
				double sigma = 0;

				if (data)
				{
					int ele = data->element(i, j, k);
					foAmp = data->data[ele][0];
					sigma = data->data[ele][1];
					free = data->getMask(i, j, k);
				}


				double diffPhwt = 0;
				double diffAmp = 0;
				
				if (diff)
				{
					diffPhwt = diff->getPhase(i, j, k);
					diffAmp = diff->getAmplitude(i, j, k);
				}


				if (foAmp != foAmp || (free == 0))
				{
					// i.e. diff of 0 when mask is free flag.
					/* No longer the case because it should be sorted before
					 * this point, and would not apply to Vagamaps */
				}

				/* MTZ file stuff */

				if (f000)
				{
					fwt = calcAmp;
					diffAmp = 0;
				}

				fdata[0] = i;
				fdata[1] = j;
				fdata[2] = k;
				fdata[3] = free;
				fdata[4] = foAmp;
				fdata[5] = sigma;
				fdata[6] = calcAmp;
				fdata[7] = fwt;
				fdata[8] = phic;
				fdata[9] = phwt;
				fdata[10] = diffAmp;
				fdata[11] = diffPhwt;

				num++;
				ccp4_lwrefl(mtzout, fdata, colout, columns, num);
			}
		}
	}

	MtzPut(mtzout, " ");
	MtzFree(mtzout);

	delete [] fdata;

}

void VagFFT::applySymmetry(bool silent)
{
	fftwf_complex *tempData;
	tempData = (fftwf_complex *)fftwf_malloc(_nn * sizeof(FFTW_DATA_TYPE));
	memset(tempData, 0, sizeof(FFTW_DATA_TYPE) * _nn);

	int count = 0;
	
	if (_spg->spg_num == 1)
	{
		return;
	}

	if (!silent)
	{
		std::cout << "applying " <<
		_spg->nsymop << " symmetry operators, space group " 
		<< _spg->symbol_xHM;
		std::cout << " (" << _spg->spg_num << ")"  << ": " << std::flush;
	}

	/* Loop through and convert data into amplitude and phase */
	for (int pre_n = 0; pre_n < _nn; pre_n++)
	{
		long n = finalIndex(pre_n);
		double xOrig = _data[n][0];
		double yOrig = _data[n][1];
		double myAmp = sqrt(xOrig * xOrig + yOrig * yOrig);
		double myPhase = atan2(yOrig, xOrig) * 180 / M_PI;
		while (myPhase >= 360) myPhase -= 360;
		while (myPhase < 0) myPhase += 360;

		_data[n][0] = myAmp;
		_data[n][1] = myPhase;
	}

	for (int k = -_nz / 2; k < _nz / 2; k++)
	{
		for (int j = -_ny / 2; j < _ny / 2; j++)
		{
			for (int i = -_nx / 2; i < _nx / 2; i++)
			{
				int abs = CSym::ccp4spg_is_sysabs(_spg, i, j, k);

				if (abs)
				{
					continue;	
				}

				long pre_index = element(i, j, k);
				long index = finalIndex(pre_index);
				/* Not misnomers: dealt with in previous paragraph */
				double myAmp = _data[index][0];
				double myPhase = _data[index][1];

				for (int l = 0; l < _spg->nsymop; l++)
				{
					float *rot = &_spg->invsymop[l].rot[0][0];

					/* rotation */
					int _h, _k, _l;
					_h = (int) rint(i*rot[0] + j*rot[3] + k*rot[6]);
					_k = (int) rint(i*rot[1] + j*rot[4] + k*rot[7]);
					_l = (int) rint(i*rot[2] + j*rot[5] + k*rot[8]);

					long sym_index = element(_h, _k, _l);
					long index = finalIndex(pre_index);
					/* translation */
					float *trn = _spg->symop[l].trn;

					double shift = (float)_h * trn[0];
					shift += (float)_k * trn[1];
					shift += (float)_l * trn[2];
					shift = fmod(shift, 1.);

					double deg = myPhase + shift * 360.;
					double newPhase = deg2rad(deg);
					
					/* add to temporary data array */
					double x = myAmp * cos(newPhase);
					double y = myAmp * sin(newPhase);

					tempData[sym_index][0] += x;
					tempData[sym_index][1] += y;
				}
			}
		}
	}

	if (!silent)
	{
		std::cout << "... done." << std::endl;
	}

	/* Loop through and convert data into amplitude and phase */
	for (int n = 0; n < _nn; n++)
	{
		long m = finalIndex(n);
		_data[m][0] = tempData[n][0];
		_data[m][1] = tempData[n][1];
	}

	free(tempData);

}

void VagFFT::multiplyAll(float value)
{
	for (long i = 0; i < nn(); i++)
	{
		int index = finalIndex(i);
		_data[index][0] *= value;
		_data[index][1] *= value;
	}
}

void VagFFT::copyFrom(VagFFTPtr other)
{
	for (long i = 0; i < nn(); i++)
	{
		int index = finalIndex(i);
		int iother = other->finalIndex(i);
		_data[index][0] = other->_data[iother][0];
		_data[index][1] = other->_data[iother][1];
	}
}

void VagFFT::copyToScratch(int scratch)
{
	for (long i = 0; i < nn(); i++)
	{
		int index = finalIndex(i);
		int sindex = scratchIndex(i, scratch);
		_data[sindex][0] = _data[index][0];
		_data[sindex][1] = _data[index][1];
	}
}

void VagFFT::addScratchBack(int scratch)
{
	for (long i = 0; i < nn(); i++)
	{
		int index = finalIndex(i);
		int sindex = scratchIndex(i, scratch);
		_data[index][0] += _data[sindex][0];
		_data[index][1] += _data[sindex][1];

	}
}
