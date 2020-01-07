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
#include <cstring>
#include <iostream>
#include <stdlib.h>

std::vector<FFTDim> VagFFT::_dimensions;

VagFFT::VagFFT(int nx, int ny, int nz, int nele)
{
	_nx = nx;
	_ny = ny;
	_nz = nz;

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
	
	return true;
}

void VagFFT::wipe()
{
	memset(_data, 0, _total);
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
	if (transform == FFTAtomsToReciprocal || transform == FFTAtomsToReal)
	{
		if (_status != FFTSeparateAtoms || _status != FFTEmpty)
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
		if (_status != FFTRealSpace || _status != FFTEmpty)
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
		if (_status != FFTReciprocalSpace || _status != FFTEmpty)
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
