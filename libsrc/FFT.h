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
#include "mat3x3.h"

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

	void prepareAtomSpace();
	void addAtom(AtomPtr atom);
	
	/** a, b, c in Angstroms; alpha, beta, gamma in degrees */
	void setUnitCell(std::vector<double> dims);
	void printSlice(int zVal = -1, double scale = 1);
private:
	long element(long x, long y, long z)
	{
		collapse(&x, &y, &z);

		return x + _nx*y + (_nx*_ny)*z;
	}

	void collapse(long *x, long *y, long *z)
	{
		while (*x < 0) *x += _nx;
		while (*x >= _nx) *x -= _nx;

		while (*y < 0) *y += _ny;
		while (*y >= _ny) *y -= _ny;

		while (*z < 0) *z += _nz;
		while (*z >= _nz) *z -= _nz;
	}

	double getAmplitude(int x, int y, int z);

	double getAmplitude(ElementPtr ele, int i, int j, int k);
	int whichColumn(ElementPtr ele);
	double cubic_interpolate(vec3 vox000, size_t im = false);
	void addInterpolatedToReal(ElementPtr ele, double sx, double sy, 
	                           double sz, double val); 
	void addExplicitAtom(AtomPtr atom);
	void addImplicitAtom(AtomPtr atom);
	/** returns true if sane */
	bool sanityCheck();

	double populateImplicit(ElementPtr ele, vec3 centre, vec3 maxVals,
	                        mat3x3 tensor, double scale, bool add);

	/** pre-loaded atom distributions converted to real space in final
	 *  column */
	void separateAtomTransform();

	void setupElements(bool wipe = true);

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
	
	mat3x3 _toReal;
	mat3x3 _toRecip;

	/* small numbers */
	mat3x3 _realBasis;

	/* big numbers */
	mat3x3 _recipBasis;
	vec3 _origin;

	bool _setMatrices;
};

#endif
