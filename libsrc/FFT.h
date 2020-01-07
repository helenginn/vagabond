// Vagabond
// Copyright (C) 2019-2020 Helen Ginn
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

#ifndef __vagabond__fft__
#define __vagabond__fft__

#include <fftw3.h>
#include "shared_ptrs.h"

typedef enum
{
	FFTEmpty,
	FFTRealSpace,
	FFTReciprocalSpace,
	FFTSeparateAtoms,
} FFTStatus;

typedef enum
{
	FFTAtomsToReal,
	FFTAtomsToReciprocal,
	FFTRealToReciprocal,
	FFTReciprocalToReal
} FFTTransform;

typedef struct
{
	int nx;
	int ny;
	int nz;
	int nele;
	fftwf_plan atom_real_to_recip;
	fftwf_plan atom_recip_to_real;
	fftwf_plan real_to_recip;
	fftwf_plan recip_to_real;
} FFTDim;

class VagFFT
{
public:
	VagFFT(int nx, int ny, int nz, int nele);
	
	void addElement(ElementPtr ele);
	
	void wipe();
	void makePlans();
	void fft(FFTTransform transform);
	void multiplyFinal(float val);
	void multiplyDotty(float val);
	
private:
	/** returns true if sane */
	bool sanityCheck();

	/** pre-loaded atom distributions converted to real space in final
	 *  column */
	void separateAtomTransform();

	std::vector<ElementPtr> _elements;
	
	FFTStatus _status;

	int _nx, _ny, _nz;
	long _nn;
	
	/** number of elements */
	int _nele;
	int _stride;

	long _total;
	fftwf_complex *_data;
	
	static std::vector<FFTDim> _dimensions;
	FFTDim *_myDims;
};

#endif
