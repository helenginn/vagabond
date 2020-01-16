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
#include "Anisotropicator.h"
#include "Model.h"
#include <cstring>
#include <float.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>

std::vector<FFTDim> VagFFT::_dimensions;

VagFFT::VagFFT(int nx, int ny, int nz, int nele)
{
	_nx = nx;
	_ny = ny;
	_nz = nz;
	_myDims = NULL;

	if (_nx % 2 == 1) _nx -= 1;
	if (_ny % 2 == 1) _ny -= 1;
	if (_nz % 2 == 1) _nz -= 1;

	_nele = nele;
	_nn = _nx * _ny * _nz;
	_status = FFTEmpty;

	_stride = 2 * (_nele) + 1;
	_total = _nn * _stride * sizeof(fftwf_complex);
	_data = (fftwf_complex *)fftwf_malloc(_total);

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
	_elements.push_back(ele);
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
		" not match" << std::endl;
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
		if (_dimensions[i].nx == _nx
		    && _dimensions[i].ny == _ny
		    && _dimensions[i].nz == _nz
		    && _dimensions[i].nele == _nele)
		{
			_myDims = &_dimensions[i];
			return;
		}
	}

	FFTDim dims;
	dims.nx = _nx; dims.ny = _ny; dims.nz = _nz; dims.nele = _nele;
	_dimensions.push_back(dims);
	_myDims = &_dimensions.back();
	
	int ns[3] = {_nz, _ny, _nx};
	unsigned fftw_flags = FFTW_MEASURE;

	fftwf_plan plan = fftwf_plan_many_dft(3, ns, _nele, _data, NULL, _stride, 
	                                      2, _data, 
	                                      NULL, _stride, 2, 
	                                      1, fftw_flags); 
	_myDims->atom_real_to_recip = plan;

	plan = fftwf_plan_many_dft(3, ns, _nele, _data, NULL, _stride, 
	                                      2, _data, 
	                                      NULL, _stride, 2, 
	                                      -1, fftw_flags); 
	_myDims->atom_recip_to_real = plan;

	int start = 2 * _nele;
	plan = fftwf_plan_many_dft(3, ns, 1, &_data[start], NULL, 
	                           _stride, 2, _data, NULL, _stride, 
	                           2, 1, fftw_flags); 

	_myDims->real_to_recip = plan;

	plan = fftwf_plan_many_dft(3, ns, 1, &_data[start], NULL, 
	                           _stride, 2, _data, NULL, _stride, 
	                           2, -1, fftw_flags); 

	_myDims->recip_to_real = plan;
}

void VagFFT::multiplyDotty(float val)
{
	for (int i = 0; i < _nn; i++)
	{
		for (int j = 0; j < _nele; j++)
		{
			long dotty_index = i * _stride + (j * 2);
			_data[dotty_index][0] *= val;
			_data[dotty_index][1] *= val;
		}
	}
}

void VagFFT::multiplyFinal(float val)
{
	for (int i = 0; i < _nn; i++)
	{
		long final_index = (i + 1) * _stride - 1;
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
					long ele_index = ele * _stride + (j * 2) + 1;
					
					ElementPtr ele = _elements[j];
					double val = 0;
					val = ele->getVoxelValue(&*ele, real.x, real.y, real.z);

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
			for (int j = 0; j < _nele; j++)
			{
				long dotty_index = i * _stride + (j * 2);

				_data[dotty_index][0] = 0; 
				_data[dotty_index][1] = 0;
			}
		}
	}

	_status = FFTSeparateAtoms;
}

void VagFFT::separateAtomTransform()
{
	fftwf_execute_dft(_myDims->atom_real_to_recip, _data, _data);
	
	for (int i = 0; i < _nn; i++)
	{
		for (int j = 0; j < _nele; j++)
		{
			long dotty_index = i * _stride + (j * 2);
			long ele_index = dotty_index + 1;

			_data[dotty_index][0] *= _data[ele_index][0];
			_data[dotty_index][1] *= _data[ele_index][1];
		}
	}

	fftwf_execute_dft(_myDims->atom_recip_to_real, _data, _data);
	multiplyDotty(1 / (double)_nn);

	for (int i = 0; i < _nn; i++)
	{
		long final_index = (i + 1) * _stride - 1;

		for (int j = 0; j < _nele; j++)
		{
			long dotty_index = i * _stride + (j * 2);

			_data[final_index][0] += _data[dotty_index][0];
			_data[dotty_index][1] += _data[dotty_index][1];
		}
	}
	
	_status = FFTRealSpace;
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
			exit(1);
		}

		separateAtomTransform();

		if (transform == FFTAtomsToReciprocal)
		{
			fftwf_execute_dft(_myDims->real_to_recip, _data, _data);
			multiplyFinal(1 / (double)_nn);
			_status = FFTReciprocalSpace;
		}
	}

	if (transform == FFTRealToReciprocal)
	{
		if ((!_status == FFTRealSpace || _status == FFTEmpty))
		{
			std::cout << "Asking for transform Real to Reciprocal, but "\
			"FFT status is not Real Space" << std::endl;
			exit(1);
		}

		fftwf_execute_dft(_myDims->real_to_recip, _data, _data);
		multiplyFinal(1 / (double)_nn);
		_status = FFTReciprocalSpace;
	}

	if (transform == FFTReciprocalToReal)
	{
		if (!(_status == FFTReciprocalSpace || _status == FFTEmpty))
		{
			std::cout << "Asking for transform Reciprocal to Real, but "\
			"FFT status is not Reciprocal Space" << std::endl;
			exit(1);
		}

		fftwf_execute_dft(_myDims->recip_to_real, _data, _data);
		multiplyFinal(1 / (double)_nn);
		_status = FFTRealSpace;
	}

}

void VagFFT::setUnitCell(std::vector<double> dims)
{
	_toReal = mat3x3_from_unit_cell(dims[0], dims[1], dims[2], dims[3],
	                                dims[4], dims[5]);
	_toRecip = mat3x3_inverse(_toReal);
	_toRecip = mat3x3_transpose(_toRecip);

	_realBasis = _toReal;
	mat3x3_scale(&_realBasis, 1 / (double)_nx, 1 / (double)_ny, 
	             1 / (double)_nz);
	_recipBasis = mat3x3_inverse(_realBasis);
	_setMatrices = true;
}

void VagFFT::addImplicitAtom(AtomPtr atom)
{
	mat3x3 tensor = atom->getModel()->getRealSpaceTensor();
	mat3x3 real_space = mat3x3_inverse(tensor);

	Anisotropicator aniso;
	aniso.setTensor(tensor);
	mat3x3 ellipsoid = aniso.basis();

	vec3 centre = atom->getAbsolutePosition();
	
	/* TODO find limiting voxels */
	const int scale = 2.5; /* no. standard devs we care about */

//	mat3x3_scale(&ellipsoid, scale, scale, scale);

	vec3 maxVals = make_vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	
	double dampen = 0;
	for (int i = 0; i < 3; i++)
	{
		vec3 axis = mat3x3_axis(ellipsoid, i);
		double length = vec3_length(axis);
		vec3_max_each(&maxVals, axis);
		vec3_mult(&axis, -1);
		vec3_max_each(&maxVals, axis);
		dampen += length / 3;
		
		/* reduction in height of gaussian with this stdev */
	}
	
	vec3_mult(&maxVals, 2.5);
	maxVals.x++; maxVals.y++; maxVals.z++;
	mat3x3_mult_vec(_recipBasis, &maxVals);
	
	vec3_subtract_from_vec3(&centre, _origin);
	mat3x3_mult_vec(_recipBasis, &centre);
	ElementPtr ele = atom->getElement();

	double total = populateImplicit(ele, centre, maxVals, tensor, 1., false);
	total = populateImplicit(ele, centre, maxVals, tensor, total, true);
}

double VagFFT::populateImplicit(ElementPtr ele, vec3 centre, vec3 maxVals,
                                mat3x3 tensor, double scale, bool add)
{
	double total = 0;
	
	for (int z = centre.z - maxVals.z; z <= centre.z + maxVals.z; z++)
	{
		for (int y = centre.y - maxVals.y; y <= centre.y + maxVals.y; y++)
		{
			for (int x = centre.x - maxVals.x; x <= centre.x + maxVals.x; x++)
			{
				vec3 pos = make_vec3(x, y, z);
				/* in voxel system */
				vec3 offset = vec3_subtract_vec3(pos, centre);
				mat3x3_mult_vec(_realBasis, &offset);

				vec3 copy = offset;
				mat3x3_mult_vec(tensor, &copy);

				copy.x *= offset.x;
				copy.y *= offset.y;
				copy.z *= offset.z;
				double multByTranspose = copy.x + copy.y + copy.z;
				double density = exp((2 * M_PI * M_PI) * -(multByTranspose));

				if (add)
				{
					addInterpolatedToReal(ele, x, y, z, density * scale);
				}
				
				total += density * scale;
			}
		}
	}

	return 1 / total;
}

int VagFFT::whichColumn(ElementPtr ele)
{
	for (int i = 0; i < _elements.size(); i++)
	{
		if (_elements[i] == ele)
		{
			return i * 2;
		}
	}

	return -1;
}

void VagFFT::addInterpolatedToReal(ElementPtr ele, double sx, double sy, 
                                   double sz, double val)
{
	long lx = (int)floor(sx);
	long ly = (int)floor(sy);
	long lz = (int)floor(sz);

	double xProps[2];
	double yProps[2];
	double zProps[2];
	double rubbish;

	xProps[1] = modf(sx + 1., &rubbish);
	yProps[1] = modf(sy + 1., &rubbish);
	zProps[1] = modf(sz + 1., &rubbish);

	xProps[0] = 1 - xProps[1];
	yProps[0] = 1 - yProps[1];
	zProps[0] = 1 - zProps[1];
	
	int column = whichColumn(ele);
	
	if (column < 0)
	{
//		std::cout << "No element!" << std::endl;
		return;
	}

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

double VagFFT::getAmplitude(ElementPtr ele, int x, int y, int z)
{
	int col = whichColumn(ele);
	long index = element(x, y, z) * _stride + col;

	double val =  (_data[index][0] * _data[index][0] + 
	               _data[index][1] * _data[index][1]);
	
	return sqrt(val);
}

double VagFFT::getAmplitude(int x, int y, int z)
{
	long index = (1 + element(x, y, z)) * _stride - 1;

	double val =  (_data[index][0] * _data[index][0] + 
	               _data[index][1] * _data[index][1]);
	
	return sqrt(val);
}

void VagFFT::printSlice(int zVal, double scale)
{
	sanityCheck();

	if (scale < 0)
	{
		scale = 1;//averageAll() * 5;
	}
	for (int j = 0; j < _ny; j++)
	{
		std::cout << "| ";
		for (int i = 0; i < _nx; i++)
		{
			std::string symbol = " ";
			double value = getAmplitude(i, j, zVal);

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
	vec3 remain = make_vec3(vox000.x - (double)((int)vox000.x),
	                        vox000.y - (double)((int)vox000.y),
	                        vox000.z - (double)((int)vox000.z));
	
	double uvw[3] = {remain.x, remain.y, remain.z};

	/* Extra refers to the additional index to be fished
	 * for 11-point interpolation. We already get 0 and 1. */
	int extra[3] = {-1, -1, -1};
	double central[3] = {0, 0, 0};
	double next[3] = {1, 1, 1};
	
	/* If uvw components are greater than 0.5, then flip them 
	 * make the extra index one ahead and reverse the order */
	for (int i = 0; i < 3; i++)
	{
		if (uvw[i] > 0.5)
		{
			extra[i] = 2;
			central[i] = 1;
			next[i] = 0;
			uvw[i] = 1 - uvw[i];
		}
	}

	long vox000x = vox000.x + central[0];
	long vox000y = vox000.y + central[1];
	long vox000z = vox000.z + central[2];
	long vox000xm = vox000.x + next[0];
	long vox000ym = vox000.y + next[1];
	long vox000zm = vox000.z + next[2];
	long vox000xn = vox000.x + extra[0];
	long vox000yn = vox000.y + extra[1];
	long vox000zn = vox000.z + extra[2];

	collapse(&vox000x, &vox000y, &vox000z);
	collapse(&vox000xm, &vox000ym, &vox000zm);
	collapse(&vox000xn, &vox000yn, &vox000zn);

	vox000y  *= _nx;
	vox000ym *= _nx;
	vox000yn *= _nx;
	vox000z  *= _nx * _ny;
	vox000zm *= _nx * _ny;
	vox000zn *= _nx * _ny;

	long int idx000 = vox000x + vox000y + vox000z;
	long int idx100 = vox000xm + vox000y + vox000z;
	long int idx010 = vox000x + vox000ym + vox000z;
	long int idx110 = vox000xm + vox000ym + vox000z;
	long int idx001 = vox000x + vox000y + vox000zm;
	long int idx101 = vox000xm + vox000y + vox000zm;
	long int idx011 = vox000x + vox000ym + vox000zm;
	long int idx111 = vox000xm + vox000ym + vox000zm;
	
	long int idxn00 = vox000xn + vox000y + vox000z;
	long int idx0n0 = vox000x + vox000yn + vox000z;
	long int idx00n = vox000x + vox000y + vox000zn;
	
	double u = uvw[0];
	double v = uvw[1];
	double w = uvw[2];
	
	double p000 = _data[(idx000 + 1) * _stride - 1][im];
	double p001 = _data[(idx001 + 1) * _stride - 1][im];
	double p010 = _data[(idx010 + 1) * _stride - 1][im];
	double p011 = _data[(idx011 + 1) * _stride - 1][im];
	double p100 = _data[(idx100 + 1) * _stride - 1][im];
	double p101 = _data[(idx101 + 1) * _stride - 1][im];
	double p110 = _data[(idx110 + 1) * _stride - 1][im];
	double p111 = _data[(idx111 + 1) * _stride - 1][im];
	
	double a = p100 - p000;
	double b = p010 - p000;
	double c = p110 - p010;
	double d = p101 - p001;
	
	double pn00 = _data[(idxn00 + 1) * _stride - 1][im];
	double p0n0 = _data[(idx0n0 + 1) * _stride - 1][im];
	double p00n = _data[(idx00n + 1) * _stride - 1][im];

	double p8value = p000+u*(a+w*(-a+d)+v*((c-a)+w*( a-c-d-p011+p111)))
	+ v*(b+w*(-p001+p011-b))+w*(-p000+p001);
	
	double mod = (p000 - 0.5 * p100 - 0.5 * pn00) * (u - u * u);
	mod += (p000 - 0.5 * p010 - 0.5 * p0n0) * (v - v * v);
	mod += (p000 - 0.5 * p001 - 0.5 * p00n) * (w - w * w);
	
	double p11value = p8value + 0.4 * mod;

	return p11value;
}

