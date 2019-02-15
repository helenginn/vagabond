// Vagabond
// Copyright (C) 2017-2018 Helen Ginn, 2011 Anton Barty
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "fftw3d.h"
#include "mat3x3.h"
#include "vec3.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "maths.h"
#include "FileReader.h"
#include "Options.h"
#include <sys/stat.h>

#include "../libccp4/cmtzlib.h"
#include "../libccp4/csymlib.h"
#include "../libccp4/ccp4_spg.h"
#include "../libccp4/ccp4_general.h"

std::deque<FourierDimension> FFT::_dimensions;

FFT::FFT()
{
	_setupBlurring = false;
	nx = 0;
	ny = 0;
	nz = 0;
	nn = 0;
	data = NULL;
	mask = NULL;
	_myDims = NULL;
	_writeToMaskZero = true;
}


FFT::FFT(long n)
{
	create(n);
}

FFT::~FFT()
{
	free(mask);
	mask = NULL;

	fftwf_free(data);
	data = NULL;
	/*
	fftwf_cleanup_threads();
	*/
}

void FFT::cleanupPlans()
{
	for (int i = 0; i < _dimensions.size(); i++)
	{
		//		fftwf_destroy_plan(_dimensions[i].plan);
		//		fftwf_destroy_plan(_dimensions[i].iplan);
	} 

	fftwf_cleanup();
}

void FFT::create(long n)
{
	create(n,n,n);
}

FFT::FFT(FFT &other)
{
	_setupBlurring = false;
	nx = other.nx;
	ny = other.ny;
	nz = other.nz;
	nn = other.nn;

	memcpy(scales, other.scales, 3 * sizeof(double));
	mask = NULL;
	data = NULL;

	if (other.data)
	{
		data = (FFTW_DATA_TYPE *)fftwf_malloc(nn * sizeof(FFTW_DATA_TYPE));
		memcpy(data, other.data, nn * sizeof(FFTW_DATA_TYPE));
	}

	if (other.mask)
	{
		mask = (int *)malloc(nn * sizeof(int));
		
		if (mask)
		{
			memcpy(mask, other.mask, nn * sizeof(int));
		}
	}

	_myDims = other._myDims;

	_basis = other._basis;
	_inverse = other._inverse;
	_writeToMaskZero = other._writeToMaskZero;
}

void FFT::copyFrom(FFTPtr other)
{
	for (int i = 0; i < nn; i++)
	{
		data[i][0] = other->data[i][0];
	}
}


void FFT::create(long nnx, long nny, long nnz)
{
	if (nnx % 2 == 1) nnx -= 1;
	if (nny % 2 == 1) nny -= 1;
	if (nnz % 2 == 1) nnz -= 1;
	nx = nnx;
	ny = nny;
	nz = nnz;
	nn = nx*ny*nz;

	if (data)
	{
		free(data);
		data = NULL;
	}

	data = (FFTW_DATA_TYPE*) fftwf_malloc(nn*sizeof(FFTW_DATA_TYPE));
	if (!data)
	{
		printf("ERROR in fftwData: Malloc failed for nn = %i\n", nn);
		exit(1);
	}

	memset(data, 0, sizeof(FFTW_DATA_TYPE) * nn);
}

double FFT::sumImag()
{
	double sum = 0;
	
	for (int i = 0; i < nn; i++)
	{
		sum += fabs(data[i][1]);
	}
	
	return sum;
}

void FFT::scaleToFFT(FFTPtr other)
{
	double ave1 = averageAll();
	double ave2 = other->averageAll();
	
	double scale = ave2 / ave1;
	
	multiplyAll(scale);
}

void FFT::setupMask()
{
	if (mask)
	{
		free(mask);
		mask = NULL;
	}

	mask = (int *) calloc(nn, sizeof(int));
}

/*
 *    Convert 3D (xyz) triple into 1D array index
 *    x is fastest axis, z is slowest axis
 */

void FFT::collapseFrac(double *xfrac, double *yfrac, double *zfrac)
{
	while (*xfrac < 0) *xfrac += 1;
	while (*xfrac >= 1) *xfrac -= 1;

	while (*yfrac < 0) *yfrac += 1;
	while (*yfrac >= 1) *yfrac -= 1;

	while (*zfrac < 0) *zfrac += 1;
	while (*zfrac >= 1) *zfrac -= 1;
}

long FFT::elementFromFrac(double xfrac, double yfrac, double zfrac)
{
	return elementFromUncorrectedFrac(xfrac, yfrac, zfrac);

	collapseFrac(&xfrac, &yfrac, &zfrac);

	double x = xfrac * nx;
	double y = yfrac * ny;
	double z = zfrac * nz;

	long index = element(x, y, z);

	return index;
}

long FFT::elementFromUncorrectedFrac(double xfrac, double yfrac, double zfrac)
{
	collapseFrac(&xfrac, &yfrac, &zfrac);

	double x = xfrac * nx;
	double y = yfrac * ny;
	double z = zfrac * nz;

	long index = element(lrint(x), lrint(y), lrint(z));

	return index;
}

void FFT::shiftToCentre()
{
	int x1, y1, z1;
	int e0, e1;
	int sx = nx / 2;
	int sy = ny / 2;
	int sz = nz / 2;

	int copyLength = sx;
	fftwf_complex *temp = (FFTW_DATA_TYPE *)fftwf_malloc(nn*sizeof(FFTW_DATA_TYPE));
	memset(temp, 0, nn*sizeof(FFTW_DATA_TYPE));
	int *tmpMask = NULL;

	if (mask)
	{
		tmpMask = (int *)calloc(nn, sizeof(int));
	}

	for (int z0 = 0; z0 < nz; z0++)
	{
		z1 = z0 + sz;
		if (z1 < 0) z1 += nz;
		if (z1 >= nz) z1 -= nz;

		for (int y0 = 0; y0 < ny; y0++)
		{
			y1 = y0 + sy;
			if (y1 < 0) y1 += ny;
			if (y1 >= ny) y1 -= ny;

			for (int x0 = 0; x0 < nx; x0 += copyLength)
			{
				x1 = x0 + sx;
				if(x1 < 0) x1 += nx;
				if(x1 >= nx) x1 -= nx;

				e0 = element(x0,y0,z0);
				e1 = element(x1,y1,z1);

				int size = sizeof(FFTW_DATA_TYPE) * copyLength;

				memcpy(&temp[e0], &data[e1], size);
				if (mask)
				{
					int maskSize = sizeof(int) * copyLength;
					memcpy(&tmpMask[e0], &mask[e1], maskSize);
				}
			}
		}
	}

	fftwf_free(data);
	data = temp;

	if (mask)
	{
		free(mask);
		mask = tmpMask;
	}
}


int FFT::setTotal(float newTotal)
{
	float total = 0;

	for (int i = 0; i < nn; i++)
	{
		total += data[i][0];
	}

	if (total <= 0) return 1;
	float scale = newTotal / total;

	for (int i = 0; i < nn; i++)
	{
		data[i][0] *= scale;
	}

	return 0;
}

void FFT::normalise()
{
	double total = averageAll();
	double mult = 1 / total;
	multiplyAll(mult);
}

void FFT::valueMinus(float value)
{
	for(long i=0; i<nn; i++)
	{
		data[i][0] = value - data[i][0];
	}
}

void FFT::cap(float value)
{
	for(long i=0; i<nn; i++)
	{
		if (data[i][0] > value)
		{
			data[i][0] = value;
		}
	}
}

void FFT::setAll(float value)
{
	for(long i=0; i<nn; i++)
	{
		data[i][0] = value;
		data[i][1] = value;
	}
}

void FFT::nxyzFromElement(long int element, long *x, long *y, long *z)
{
	*x = element % nx;
	element -= *x;
	element /= nx;

	*y = element % ny;
	element -= *y;
	element /= ny;

	*z = element;
}

vec3 FFT::fracFromElement(long int element)
{
	long x = element % nx;
	element -= x;
	element /= nx;

	long y = element % ny;
	element -= y;
	element /= ny;

	long z = element;

	double xfrac = (double)x / (double)nx;
	double yfrac = (double)y / (double)ny;
	double zfrac = (double)z / (double)nz;

	return make_vec3(xfrac, yfrac, zfrac);
}

// in degrees.
double FFT::getPhase(long x, long y, long z)
{
	long index = element(x, y, z);

	double degrees = atan2(data[index][1], data[index][0]) * 180 / M_PI;

	while (degrees >= 360) degrees -= 360;

	while (degrees < 0) degrees += 360;

	return degrees;
}

double FFT::getIntensity(long index)
{
	return (data[index][0] * data[index][0] + data[index][1] * data[index][1]);

}

double FFT::getIntensity(long x, long y, long z)
{
	long index = element(x, y, z);

	return (data[index][0] * data[index][0] + data[index][1] * data[index][1]);
}

double FFT::getImaginary(long x, long y, long z)
{
	long index = element(x, y, z);

	return data[index][1];
}

void FFT::setReal(double xfrac, double yfrac, double zfrac, double real)
{
	long index = elementFromFrac(xfrac, yfrac, zfrac);

	data[index][0] = real;
	data[index][1] = 0;
}

void FFT::addToReal(double xfrac, double yfrac, double zfrac, double real)
{
	long index = elementFromFrac(xfrac, yfrac, zfrac);

	data[index][0] += real;
}

#define START_LOOP -1
#define END_LOOP 2

void FFT::setupBlurring()
{
	_blurAmounts.clear();

	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	double bfac = crystal->getRealBFactor();
	bfac /= 8 * M_PI * M_PI;

	for (int i = START_LOOP; i < END_LOOP; i++)
	{
		for (int j = START_LOOP; j < END_LOOP; j++)
		{
			for (int k = START_LOOP; k < END_LOOP; k++)
			{
				vec3 shift = make_vec3(i, j, k);
				mat3x3_mult_vec(_basis, &shift);
				double movement = vec3_sqlength(shift);

				float factor = normal_distribution(movement, bfac);

				_blurAmounts.push_back(factor);
			}
		}
	}		

	_setupBlurring = true;
}

void FFT::addInterpolatedToReal(double sx, double sy, double sz, double val)
{
	long lx = (int)floor(sx);
	long ly = (int)floor(sy);
	long lz = (int)floor(sz);

	double xProps[2];
	double yProps[2];
	double zProps[2];

	xProps[1] = fmod(sx + 1, 1.);
	yProps[1] = fmod(sy + 1, 1.);
	zProps[1] = fmod(sz + 1, 1.);

	xProps[0] = 1 - xProps[1];
	yProps[0] = 1 - yProps[1];
	zProps[0] = 1 - zProps[1];

	for (int p = 0; p < 2; p++)
	{
		for (int q = 0; q < 2; q++)
		{
			for (int r = 0; r < 2; r++)
			{
				int sx1 = lx + p;
				int sy1 = ly + q;
				int sz1 = lz + r;

				long index = element(sx1, sy1, sz1);
				double prop = xProps[p] * yProps[q] * zProps[r];
				data[index][0] += prop * val;
			}	
		}
	}
}

void FFT::addInterpolatedToFrac(double fx, double fy, double fz, double val)
{
	collapseFrac(&fx, &fy, &fz);

	double sx = fx * nx;
	double sy = fy * ny;
	double sz = fz * nz;

	addInterpolatedToReal(sx, sy, sz, val);
}

void FFT::addBlurredToReal(double xfrac, double yfrac, double zfrac, double real)
{
	if (!_setupBlurring)
	{
		setupBlurring();
	}

	collapseFrac(&xfrac, &yfrac, &zfrac);

	double x = xfrac * nx;
	double y = yfrac * ny;
	double z = zfrac * nz;

	int count = 0;
	for (int k = START_LOOP; k < END_LOOP; k++)
	{
		double sz = z + (double)k;
		for (int j = START_LOOP; j < END_LOOP; j++)
		{
			double sy = y + (double)j;
			for (int i = START_LOOP; i < END_LOOP; i++)
			{
				double sx = x + (double)i;
				long lx = (int)floor(sx);
				long ly = (int)floor(sy);
				long lz = (int)floor(sz);

				long index = element(lx, ly, lz);

				if (!_writeToMaskZero && mask[index] == 0)
				{
					continue;
				}

				float factor = _blurAmounts[count];
				count++;

				addInterpolatedToReal(sx, sy, sz, real * factor);
			}
		}
	}
}

void FFT::multiplyAll(float value)
{
	for(long i = 0; i < nn; i++)
	{
		data[i][0] *= value;
		data[i][1] *= value;
	}
}

void FFT::takePlansFrom(FFTPtr other)
{
	_myDims = other->_myDims;
}

void FFT::createFFTWplan(int nthreads, unsigned fftw_flags)
{
	nthreads = 1;
	char    wisdomFile[128];
	wisdomFile[0] = 0;

	for (int i = 0; i < _dimensions.size(); i++)
	{
		if (_dimensions[i].nx == nx
		    && _dimensions[i].ny == ny
		    && _dimensions[i].nz == nz)
		{
			_myDims = &_dimensions[i];
			return;
		}
	}


	/*
	 *    Sanity check
	 */
	if (nx<=0 || ny<=0 || nz<=0)
	{
		printf("Illogical FFT dimensions: %li x %li x %li\n", nx,ny,nz);
		exit(1);
	}

	if (sizeof(FFTW_DATA_TYPE) == 2*sizeof(float))
	{
		strcat(wisdomFile,".fftw3f_wisdom");
	}
	else if (sizeof(FFTW_DATA_TYPE) == 2*sizeof(double))
	{
		strcat(wisdomFile,".fftw3_wisdom");
	}

	if (file_exists(wisdomFile))
	{
		fftwf_import_wisdom_from_filename(wisdomFile);
	}

	FourierDimension dims;
	dims.nx = nx; dims.ny = ny; dims.nz = nz;
	_dimensions.push_back(dims);
	_myDims = &_dimensions[_dimensions.size() - 1];

	/* Generate FFTW plans */
	_myDims->plan = fftwf_plan_dft_3d((int)nz, (int)ny, (int)nx, 
	                                  data, data, 1, fftw_flags);
	_myDims->iplan = fftwf_plan_dft_3d((int)nz, (int)ny, (int)nx, 
	                                   data, data, -1, fftw_flags);

	/*  Export wisdom to file */
	int success = fftwf_export_wisdom_to_filename(wisdomFile);

	if (!success)
	{
		printf("\t\tCould not export wisdom to %s\n",wisdomFile);
	}

	/* Grab it back */
	fftwf_import_wisdom_from_filename(wisdomFile);

	/* Use it */
	fft(1);
	fft(-1);

	/* Reinitialise */
	setAll(0);
}



/*
 *    Do the FFT
 */
void FFT::fft(int direction)
{
	if (direction == 1)
	{
		fftwf_execute_dft(_myDims->plan, data, data);
	}
	else if (direction == -1)
	{
		fftwf_execute_dft(_myDims->iplan, data, data);
	}
	else
	{
		printf("Error: Undefined FFTW direction\n");
		exit(1);
	}

}


void FFT::setBasis(mat3x3 mat, double sampleScale)
{
	_basis = mat;
	mat3x3_scale(&_basis, sampleScale, sampleScale, sampleScale);
	scales[0] = mat3x3_length(_basis, 0);
	scales[1] = mat3x3_length(_basis, 1);
	scales[2] = mat3x3_length(_basis, 2);

	_inverse = mat3x3_inverse(_basis);
}

void FFT::invertScale()
{
	setBasis(_inverse, 1);
}

/* 11-point interpolation - attempted transcription from Dave's function
 * from GAP */
double FFT::cubic_interpolate(vec3 vox000, size_t im)
{
	/* vox000 has integer and real components */
	
	/* Pick out just the real components */
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

	vox000y  *= nx;
	vox000ym *= nx;
	vox000yn *= nx;
	vox000z  *= nx * ny;
	vox000zm *= nx * ny;
	vox000zn *= nx * ny;

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
	
	double p000 = data[idx000][im];
	double p001 = data[idx001][im];
	double p010 = data[idx010][im];
	double p011 = data[idx011][im];
	double p100 = data[idx100][im];
	double p101 = data[idx101][im];
	double p110 = data[idx110][im];
	double p111 = data[idx111][im];
	
	double a = p100 - p000;
	double b = p010 - p000;
	double c = p110 - p010;
	double d = p101 - p001;
	
	double pn00 = data[idxn00][im];
	double p0n0 = data[idx0n0][im];
	double p00n = data[idx00n][im];

	double p8value = p000+u*(a+w*(-a+d)+v*((c-a)+w*( a-c-d-p011+p111)))
	+ v*(b+w*(-p001+p011-b))+w*(-p000+p001);
	
	double mod = (p000 - 0.5 * p100 - 0.5 * pn00) * (u - u * u);
	mod += (p000 - 0.5 * p010 - 0.5 * p0n0) * (v - v * v);
	mod += (p000 - 0.5 * p001 - 0.5 * p00n) * (w - w * w);
	
	double p11value = p8value + 0.4 * mod;

	return p11value;
}

double FFT::score(FFTPtr fftCrystal, FFTPtr fftThing, vec3 pos,
                  std::vector<CoordVal> *vals, MapScoreType mapScore)
{
	return operation(fftCrystal, fftThing, pos, mapScore, vals);
}

void FFT::findLimitingValues(double xMin, double xMax, double yMin,
                             double yMax, double zMin, double zMax,
                             vec3 *minVals, vec3 *maxVals)
{
	mat3x3 toCrystBasis = getBasisInverse();
	*minVals = make_vec3(FLT_MAX, FLT_MAX, FLT_MAX);;
	*maxVals = make_vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);;

	for (double k = zMin; k <= zMax + 0.1; k += (zMax - zMin))
	{
		for (double j = yMin; j <= yMax + 0.1; j += (zMax - zMin))
		{
			for (double i = xMin; i <= xMax + 0.1; i += (zMax - zMin))
			{
				vec3 angCorner = make_vec3(i, j, k);
				vec3 toCrystal = mat3x3_mult_vec(toCrystBasis, angCorner);
				toCrystal.x = lrint(toCrystal.x);
				toCrystal.y = lrint(toCrystal.y);
				toCrystal.z = lrint(toCrystal.z);

				vec3_min_each(minVals, toCrystal);
				vec3_max_each(maxVals, toCrystal);
			}
		}
	}
}

typedef struct
{
	long element;
	double value;
} IndexValue;

void FFT::blurRealToImaginary(int x, int y, int z, mat3x3 tensor)
{
	long ele = element(x, y, z);
	double val = data[ele][0];
	
	double longest = std::max(tensor.vals[0], 
	                          std::max(tensor.vals[4], tensor.vals[8]));

	vec3 mins = empty_vec3();
	vec3 maxs = empty_vec3();
	findLimitingValues(-longest, longest, -longest, longest, 
	                   -longest, longest, &mins, &maxs);
	
	mat3x3 basis = getBasis();
	
	std::vector<IndexValue> additions;
	double count = 0;

	for (int k = mins.z; k < maxs.z; k++)
	{
		for (int j = mins.y; j < maxs.y; j++)
		{
			for (int i = mins.x; i < maxs.x; i++)
			{
				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(basis, &ijk);
				mat3x3_mult_vec(tensor, &ijk);

				ijk.x *= i; ijk.y *= j; ijk.z *= k;
				double dist = ijk.x + ijk.y + ijk.z;
				double aniso = exp(2 * M_PI * M_PI * -(dist));
				
				long blele = element(x + i, y + j, z + k);
				
				IndexValue pair;
				pair.element = blele;
				pair.value = aniso * val;
				additions.push_back(pair);
				
				count += aniso * val;
				
			}
		}
	}
	
	double ratio = 1 / count;
	
	for (int i = 0; i < additions.size(); i++)
	{
		long ele = additions[i].element;
		double add = ratio * additions[i].value;
		data[ele][1] += add;
	}
}

void FFT::shrink(double radius)
{
	vec3 mins = make_vec3(0, 0, 0);
	vec3 maxs = make_vec3(0, 0, 0);
	mat3x3 basis = getBasis();
	findLimitingValues(0, radius, 0, radius, 0, radius,
					   &mins, &maxs);
					
	int count = 0;
	int total = 0;
	int solv = 0;

	for (int z = 0; z < nz; z++)
	{
		for (int y = 0; y < ny; y++)
		{
			for (int x = 0; x < nx; x++)
			{
				long raw = element(x, y, z);
				/* We only want to modify protein to become
				 * more like solvent */
				if (data[raw][0] == 1)
				{
					solv++;
					data[raw][1] = 1;
					continue;
				}

				/* Default is to be protein */
				data[raw][1] = 0;
				bool done = false;

				for (int k = 0; k < maxs.z && !done; k++)
				{
					for (int j = 0; j < maxs.y && !done; j++)
					{
						for (int i = 0; i < maxs.x && !done; i++)
						{
							vec3 ijk = make_vec3(i, j, k);
							mat3x3_mult_vec(basis, &ijk);

							/* Doesn't matter if it goes over radial
							 * boundary */
							if (vec3_sqlength(ijk) > radius * radius)
							{
								continue;
							}

							long index = element(i + x, j + y, k + z);
							
							if (data[index][0] == 1)
							{
								count++;
								done = true;
								data[raw][1] = 1;
								break;
							}
						}
					}
				}
				
				total++;
			}
		}
	}
	
	for (int i = 0; i < nn; i++)
	{
		data[i][0] = data[i][1];
		data[i][1] = 0;
	}
}

void FFT::addToValueAroundPoint(vec3 pos, double radius, double value)
{
	/* Determine square bounding box for radius in Ang. */
	vec3 minRadius = make_vec3(0, 0, 0);
	vec3 maxRadius = make_vec3(0, 0, 0);
	
	collapseFrac(&pos.x, &pos.y, &pos.z);
	pos.x *= nx;
	pos.y *= ny;
	pos.z *= nz;

	findLimitingValues(-radius, radius, -radius, radius, -radius, radius,
	                   &minRadius, &maxRadius);
	mat3x3 basis = getBasis();

	for (double k = minRadius.z; k < maxRadius.z; k++)
	{
		for (double j = minRadius.y; j < maxRadius.y; j++)
		{
			for (double i = minRadius.x; i < maxRadius.x; i++)
			{
				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(basis, &ijk);
				
				/* Doesn't matter if it goes over radial
				 * boundary */
				if (vec3_sqlength(ijk) > radius * radius)
				{
					continue;
				}

				double x = lrint(i + pos.x);
				double y = lrint(j + pos.y);
				double z = lrint(k + pos.z);
				
				long index = element(x, y, z);
				data[index][0] += value;
				
				if (data[index][0] < 0)
				{
					data[index][0] = 0;
				}
			}
		}
	}
}


/*  For multiplying point-wise
 *  No assumption that interpolation is not needed.
 */
double FFT::operation(FFTPtr fftEdit, FFTPtr fftConst, vec3 add,
                      MapScoreType mapScoreType, std::vector<CoordVal> *vals,
                      bool sameScale, bool interp)
{
	/* I rarely comment something so heavily but I will get confused if
	 * I don't, this time, as I can't soak the protocol into the variable
	 * names. Bear in mind the three coordinate systems:
	 * (a) Angstroms
	 * (b) Crystal voxels
	 * (c) Atom voxels */

	FFT *fftCrystal = &*fftEdit;
	FFT *fftAtom = &*fftConst;
	double volume = 1;

	/* Bring the fractional coordinate of the atom into range 0 < frac <= 1 */
	if (mapScoreType != MapScoreAddNoWrap)
	{
		FFT::collapseFrac(&add.x, &add.y, &add.z);
	}

	/* Multiply by the relative dimensions of the crystal */
	double multX = add.x * fftCrystal->nx;
	double multY = add.y * fftCrystal->ny;
	double multZ = add.z * fftCrystal->nz;

	/* Get the remainder after subtracting a number of whole voxels */
	vec3 atomOffset; // in crystal coordinates actually, converted later.
	atomOffset.x = fmod(multX, 1.);
	atomOffset.y = fmod(multY, 1.);
	atomOffset.z = fmod(multZ, 1.);

	/* Store the non-remainder whole voxel values for way later. */
	vec3 atomWholeCoords = make_vec3((int)multX, (int)multY, (int)multZ);

	/* Prepare a matrix to convert crystal voxels into atomic voxels */
	mat3x3 crystal2AtomVox = mat3x3_mult_mat3x3(fftAtom->getBasisInverse(),
	                                            fftCrystal->getBasis());

	/* Prepare a matrix to convert atomic voxels into crystal voxels */
	mat3x3 atomVox2Crystal = mat3x3_mult_mat3x3(fftCrystal->getBasisInverse(),
	                                            fftAtom->getBasis());

	/* Apply this offset and reverse it. This small offset must be added
	 * to all future atomic coordinates prior to interpolation. This
	 * is therefore now in atom voxels.*/
	mat3x3_mult_vec(crystal2AtomVox, &atomOffset);
	vec3_mult(&atomOffset, -1);

	fftAtom->shiftToCentre();
	
	mat3x3 crystBasis = fftCrystal->getBasis();
	mat3x3 atomBasis = fftAtom->getBasis();
	double vol_corr = mat3x3_volume(crystBasis) / mat3x3_volume(atomBasis);

	/* There will be an additional shift having moved the atom by
	 * half the dimension length which needs to be taken into account, 
	 * unfortunately. */
	vec3 atomShift = make_vec3((double)(-fftAtom->nx) * 0.5,
	                           (double)(-fftAtom->ny) * 0.5,
	                           (double)(-fftAtom->nz) * 0.5);
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
	vec3 cornerCrystal = vec3_add_vec3(wholeShiftOnly, atomWholeCoords);

	/* Fractional offset in atomic coordinates, for each atom as a
	 * 	fraction of the crystal voxel. */
	atomOffset = vec3_add_vec3(atomOffset, shiftRemainder);
	vec3 crystOffset = mat3x3_mult_vec(atomVox2Crystal, atomOffset);

	/* We loop around these crystal voxel limits now (ss -> ms -> fs).
	 * We also break the loop if it exceeds the limits of our atom voxels
	 * during the loop itself. */

	int count = 0;

	/* Determine bounding box - 9th Dec 2017 */
	vec3 minAtom = make_vec3(0, 0, 0);
	vec3 maxAtom = make_vec3(0, 0, 0);

	for (int k = 0; k <= fftAtom->nz; k += fftAtom->nz)
	{
		for (int j = 0; j <= fftAtom->ny; j += fftAtom->ny)
		{
			for (int i = 0; i <= fftAtom->nx; i += fftAtom->nx)
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

	/* Set all the voxels to zero if we are going to copy across info.
	 * We do not want any contamination with original Fc. */
	if (mapScoreType == MapScoreTypeCopyToSmaller)
	{
		fftAtom->setAll(0);
	}
	
	double step = 1;

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

					if (atomPos.x > fftAtom->nx || atomPos.y > fftAtom->ny
					    || atomPos.z > fftAtom->nz)
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
				vec3_add_to_vec3(&atomPos, atomOffset);

				/* If this value is within floating point error, stop now. */
				if (fftAtom->getReal(atomPos.x, atomPos.y, atomPos.z) <= 10e-6
				    && mapScoreType != MapScoreTypeCopyToSmaller)
				{
					continue;
				}

				/* Find the interpolated value which atomPos falls on */
				double atomReal = 0;

				if (atomPos.x < 0) atomPos.x += fftAtom->nx;
				if (atomPos.y < 0) atomPos.y += fftAtom->ny;
				if (atomPos.z < 0) atomPos.z += fftAtom->nz;

				if (interp)
				{
					atomReal = fftAtom->cubic_interpolate(atomPos, 0);
				}
				else
				{
					atomReal = fftAtom->getReal(lrint(atomPos.x),
                                                lrint(atomPos.y), 
                                                lrint(atomPos.z));
				}

				/* We add the crystal offset so we don't end up with thousands
				 * of atoms at the very centre of our map */
				vec3 cVox = vec3_add_vec3(crystalPos, cornerCrystal);

				if (mapScoreType == MapScoreAddNoWrap)
				{
					if (cVox.x < -fftCrystal->nx / 2 || 
					    cVox.y < -fftCrystal->ny / 2 ||
					    cVox.z < -fftCrystal->nz / 2)
					{
						continue;
					}					

					if (cVox.x > fftCrystal->nx / 2 ||
					    cVox.y > fftCrystal->ny / 2 ||
					    cVox.z > fftCrystal->nz / 2)
					{
						continue;
					}
				}

				if (cVox.x < 0) cVox.x += fftCrystal->nx;
				if (cVox.y < 0) cVox.y += fftCrystal->ny;
				if (cVox.z < 0) cVox.z += fftCrystal->nz;


				/* Get the index of this final crystal voxel. */
				long cIndex = fftCrystal->element(cVox.x + 0.5,
				                                  cVox.y + 0.5,
				                                  cVox.z + 0.5);

				count++;

				if (mapScoreType == MapScoreTypeCorrel)
				{
					if ((!fftCrystal->_writeToMaskZero &&
					     fftCrystal->getMask(cIndex) == 0))
					{
						continue;
					}

					// We do NOT need to interpolate //
					double realCryst = fftCrystal->getReal(cIndex);

					if (vals)
					{
						CoordVal val;
						val.fo = realCryst;
						val.fc = atomReal;
						#ifdef COORDVAL_FULL
						vec3 frac = fftCrystal->fracFromElement(cIndex);
						val.pos = frac;
						#endif

						vals->push_back(val);
					}
				}
				else if (mapScoreType == MapScoreTypeCopyToSmaller)
				{
					double realCryst = fftCrystal->getReal(cIndex);
					int ele = fftAtom->element(atomPos);

					fftAtom->setElement(ele, realCryst, 0);
				}
				else if (mapScoreType == MapScoreTypeNone ||
				         mapScoreType == MapScoreAddNoWrap)
				{
					/* Add the density to the real value of the crystal
					 * voxel. */

					if (fftCrystal->_writeToMaskZero ||
					    fftCrystal->getMask(cIndex) != 0)
					{
						fftCrystal->data[cIndex][0] += atomReal * vol_corr;
					}
				}
			}
		}
	}
	
	
	if (mapScoreType == MapScoreTypeCopyToSmaller)
	{
		fftAtom->shiftToCentre();
	}

	return 0;
}

/*  For multiplying point-wise
 *  Assumes no differences in FFT scales.
 */
void FFT::multiply(FFTPtr fftEdit, FFTPtr fftConst)
{
	FFT *fftSmall = &*fftConst;
	FFT *fftBig = &*fftEdit;

	for (long int i = 0; i < fftSmall->nn; i++)
	{
		float real = fftBig->data[i][0] * fftSmall->data[i][0]
		- fftBig->data[i][1] * fftSmall->data[i][1];
		fftBig->data[i][0] = real;

		float imag = fftBig->data[i][0] * fftSmall->data[i][1]
		+ fftBig->data[i][1] * fftSmall->data[i][0];
		fftBig->data[i][1] = imag;
	}
}

/*  For multiplying point-wise
 *  Assumes identical nx/ny/nz/scales.
 */
void FFT::addSimple(FFTPtr fftEdit, FFTPtr fftConst)
{
	for (long int i = 0; i < fftEdit->nn; i++)
	{
		if (!fftEdit->_writeToMaskZero && fftEdit->mask[i] == 0)
		{
			continue;
		}

		fftEdit->data[i][0] += fftConst->data[i][0];
	}
}

void FFT::printCinema()
{
	for (int i = 0; i < nz; i++)
	{
		printSlice(i);
		sleep(1);
	}
}

void FFT::printSlice(int zVal)
{
	for (int j = 0; j < ny; j++)
	{
		std::cout << "| ";
		for (int i = 0; i < nx; i++)
		{
			std::string symbol = " ";
			double value = sqrt(getIntensity(i, j, zVal));

			if (value > 0.01) symbol = ".";
			if (value > 0.02) symbol = ":";
			if (value > 0.04) symbol = "\"";
			if (value > 0.08) symbol = "*";
			if (value > 0.16) symbol = "x";
			if (value > 0.32) symbol = "H";
			if (value > 0.64) symbol = "#";
			if (value > 1.00) symbol = "@";

			std::cout << symbol;
		}

		std::cout << " |" << std::endl;
	}
	std::cout << std::endl;
}

vec3 FFT::collapseToRealASU(vec3 frac, CSym::CCP4SPG *spaceGroup, 
                            int *flipped)
{
	for (int l = 0; l < spaceGroup->nsymop; l++)
	{
		if (flipped && *flipped >= 0)
		{
			l = *flipped;
		}
	
		vec3 ijk = frac;
		float *trn = spaceGroup->symop[l].trn;
		float *rot = &spaceGroup->symop[l].rot[0][0];

		ijk.x = (int) rint(frac.x*rot[0] + frac.y*rot[3] + frac.z*rot[6]);
		ijk.y = (int) rint(frac.x*rot[1] + frac.y*rot[4] + frac.z*rot[7]);
		ijk.z = (int) rint(frac.x*rot[2] + frac.y*rot[5] + frac.z*rot[8]);

		ijk.x += trn[0];
		ijk.y += trn[1];
		ijk.z += trn[2];

		if (ijk.x < 0 || ijk.y < 0 || ijk.z < 0)
		{
			continue;
		}

		float *lim = spaceGroup->mapasu_zero;

		if (ijk.x > lim[0] || ijk.y > lim[1] || ijk.z > lim[2])
		{
			continue;
		}

		return ijk;
	}

	return frac;
}

void FFT::applySymmetry(CSym::CCP4SPG *spaceGroup, bool silent)
{
	fftwf_complex *tempData;
	tempData = (fftwf_complex *)fftwf_malloc(nn * sizeof(FFTW_DATA_TYPE));
	memset(tempData, 0, sizeof(FFTW_DATA_TYPE) * nn);

	int count = 0;

	//	spaceGroup = CSym::ccp4spg_load_by_ccp4_num(155);
	
	if (!silent)
	{
		std::cout << "applying symmetry, space group " << spaceGroup->symbol_xHM;
		std::cout << " (" << spaceGroup->spg_num << ")"  << ": " << std::flush;
	}

	/* Loop through and convert data into amplitude and phase */
	for (int n = 0; n < nn; n++)
	{
		double xOrig = data[n][0];
		double yOrig = data[n][1];
		double myAmp = sqrt(xOrig * xOrig + yOrig * yOrig);
		double myPhase = atan2(yOrig, xOrig) * 180 / M_PI;
		while (myPhase >= 360) myPhase-= 360;
		while (myPhase < 0) myPhase += 360;

		data[n][0] = myAmp;
		data[n][1] = myPhase;
	}

	for (int k = -nz / 2; k < nz / 2; k++)
	{
		for (int j = -ny / 2; j < ny / 2; j++)
		{
			for (int i = -nx / 2; i < nx / 2; i++)
			{
				int abs = CSym::ccp4spg_is_sysabs(spaceGroup, i, j, k);

				if (abs)
				{
					continue;	
				}

				long index = element(i, j, k);
				/* Not misnomers: dealt with in previous paragraph */
				double myAmp = data[index][0];
				double myPhase = data[index][1];

				for (int l = 0; l < spaceGroup->nsymop; l++)
				{
					float *rot = &spaceGroup->invsymop[l].rot[0][0];

					/* rotation */
					int _h, _k, _l;
					_h = (int) rint(i*rot[0] + j*rot[3] + k*rot[6]);
					_k = (int) rint(i*rot[1] + j*rot[4] + k*rot[7]);
					_l = (int) rint(i*rot[2] + j*rot[5] + k*rot[8]);

					long sym_index = element(_h, _k, _l);
					/* translation */
					float *trn = spaceGroup->symop[l].trn;

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

	memcpy(data, tempData, sizeof(FFTW_DATA_TYPE) * nn);
	free(tempData);
}

void FFT::writeReciprocalToFile(std::string filename, double maxResolution,
                                CSym::CCP4SPG *mtzspg, std::vector<double> unitCell,
                                mat3x3 real2Frac, FFTPtr data,
	                           std::vector<double> bins, 
	                           std::vector<double> ampAves) 
{
	double nLimit[3];
	nLimit[0] = nx;
	nLimit[1] = ny;
	nLimit[2] = nz;
	
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

	int columns = 10;

	float cell[6], wavelength;
	float *fdata = new float[columns];

	/* variables for symmetry */
	float rsm[192][4][4];
	char ltypex[2];

	/* variables for MTZ data structure */
	CMtz::MTZ *mtzout;
	CMtz::MTZXTAL *xtal;
	CMtz::MTZSET *set;
	CMtz::MTZCOL *colout[11];

	if (unitCell.size() < 6)
	{
		unitCell.resize(6);
		unit_cell_from_mat3x3(_basis, &unitCell[0]);
		unitCell[0] *= nx;
		unitCell[1] *= ny;
		unitCell[2] *= nz;

		mtzspg = CSym::ccp4spg_load_by_ccp4_num(1);
		real2Frac = mat3x3_from_unit_cell(&unitCell[0]);
		real2Frac = mat3x3_inverse(real2Frac);
		
	}

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
	for (int i = 0; i < mtzspg->nsymop; ++i)
	CCP4::rotandtrn_to_mat4(rsm[i], mtzspg->symop[i]);
	strncpy(ltypex, mtzspg->symbol_old, 1);
	ccp4_lwsymm(mtzout, mtzspg->nsymop, mtzspg->nsymop_prim, rsm, ltypex,
	            mtzspg->spg_ccp4_num, mtzspg->symbol_old, mtzspg->point_group);

	// then add xtals, datasets, cols
	xtal = MtzAddXtal(mtzout, "vagabond_crystal", "vagabond_project", cell);
	set = MtzAddDataset(mtzout, xtal, "Dataset", wavelength);
	colout[0] = MtzAddColumn(mtzout, set, "H", "H");
	colout[1] = MtzAddColumn(mtzout, set, "K", "H");
	colout[2] = MtzAddColumn(mtzout, set, "L", "H");
	colout[3] = MtzAddColumn(mtzout, set, "FREE", "I");
	colout[4] = MtzAddColumn(mtzout, set, "FP", "F");
	colout[5] = MtzAddColumn(mtzout, set, "FC", "F");
	colout[6] = MtzAddColumn(mtzout, set, "FWT", "F");
	colout[7] = MtzAddColumn(mtzout, set, "PHWT", "P");
	colout[8] = MtzAddColumn(mtzout, set, "DELFWT", "F");
	colout[9] = MtzAddColumn(mtzout, set, "PHDELWT", "P");

	int num = 0;

	/* symmetry issues */
	for (int k = -nLimit[2]; k < nLimit[2]; k++)
	{
		for (int j = -nLimit[1]; j < nLimit[1]; j++)
		{
			for (int i = -nLimit[0]; i < nLimit[0]; i++)
			{
				bool asu = CSym::ccp4spg_is_in_asu(mtzspg, i, j, k);
				bool f000 = (i == 0 && j == 0 && k == 0);

				if (!asu || f000)
				{
					continue;
				}

				vec3 pos = make_vec3(i, j, k);
				mat3x3_mult_vec(real2Frac, &pos);

				if (vec3_length(pos) > dStar)
				{
					continue;
				}

				double phase = getPhase(i, j, k);

				double intensity = getIntensity(i, j, k);
				double calcAmp = sqrt(intensity);

				int free = 1;
				double foInt = intensity;

				if (data)
				{
					foInt = data->getIntensity(i, j, k);
					free = data->getMask(i, j, k);
				}

				double foAmp = sqrt(foInt);
				double fofofc = 2 * foAmp - calcAmp;
				double fofc = foAmp - calcAmp;
				
				if (free == 0)
				{
					// Substitute Fcalc when free = 0 so as not
					// to pollute with data.
					fofofc = calcAmp;
				}

				if (foAmp != foAmp || (free == 0))
				{
					// i.e. diff of 0 when mask is free flag.
					fofc = 0;
				}

				/* MTZ file stuff */

				if (f000)
				{
					fofofc = calcAmp;
					fofc = 0;
				}

				fdata[0] = i;
				fdata[1] = j;
				fdata[2] = k;
				fdata[3] = free;
				fdata[4] = foAmp;
				fdata[5] = calcAmp;
				fdata[6] = fofofc;
				fdata[7] = phase;
				fdata[8] = fofc;
				fdata[9] = phase;

				num++;
				ccp4_lwrefl(mtzout, fdata, colout, columns, num);
			}
		}
	}

	MtzPut(mtzout, " ");
	MtzFree(mtzout);

	delete [] fdata;
}

vec3 getSymRelatedPosition(CrystalPtr crystal, vec3 pos, int i)
{
	CSym::CCP4SPG *spg = crystal->getSpaceGroup();
	mat3x3 real2Frac = crystal->getReal2Frac();
	mat3x3_mult_vec(real2Frac, &pos);

	float *rot = &spg->symop[i].rot[0][0];
	float *trn = spg->symop[i].trn;

	vec3 mod = empty_vec3();
	mod.x = pos.x * rot[0] + pos.y * rot[1] + pos.z * rot[2];
	mod.y = pos.x * rot[3] + pos.y * rot[4] + pos.z * rot[5];
	mod.z = pos.x * rot[6] + pos.y * rot[7] + pos.z * rot[8];
	mod.x += trn[0];
	mod.y += trn[1];
	mod.z += trn[2];
	
	mat3x3 frac2Real = crystal->getHKL2Frac();
	mat3x3_mult_vec(frac2Real, &mod);
	
	return mod;
}


vec3 FFT::getPositionInAsu(vec3 vec)
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	CSym::CCP4SPG *spg = crystal->getSpaceGroup();
	mat3x3 real2Frac = crystal->getReal2Frac();
	
	for (int i = 0; i < spg->nsymop; i++)
	{
		vec3 pos = getSymRelatedPosition(crystal, vec, i);
		vec3 tmp = pos;
		mat3x3_mult_vec(real2Frac, &tmp);
		FFT::collapseFrac(&tmp.x, &tmp.y, &tmp.z);
		
		if (tmp.x < spg->mapasu_zero[0] &&
		    tmp.y < spg->mapasu_zero[1] &&
		    tmp.z < spg->mapasu_zero[2]) 
		{
			return pos;
		}
	}
	
	return empty_vec3();
}

