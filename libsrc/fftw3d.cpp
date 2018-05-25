/*
*  fftw.cpp
*    A simple object wrapper for 3D FFTs
*
*  Created by Anton Barty on 26/07/11.
*  Copyright 2011 Anton Barty, 2017 Helen Ginn. All rights reserved.
*
*/

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
#include <sys/stat.h>

#include "../libccp4/cmtzlib.h"
#include "../libccp4/csymlib.h"
#include "../libccp4/ccp4_spg.h"
#include "../libccp4/ccp4_general.h"

std::deque<FourierDimension> FFT::_dimensions;

inline void fftwf_add(fftwf_complex comp1, fftwf_complex comp2, float *result)
{
	result[0] = comp1[0];
	result[1] = comp1[1];
	result[0] += comp2[0];
	result[1] += comp2[1];
}


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
	fftwf_cleanup_threads();
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
		memcpy(mask, other.mask, nn * sizeof(int));
	}

	_myDims = other._myDims;

	_basis = other._basis;
	_inverse = other._inverse;
	_writeToMaskZero = other._writeToMaskZero;
}

void FFT::copyFrom(FFTPtr other)
{
	memcpy(data, other->data, nn * sizeof(FFTW_DATA_TYPE));
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

void FFT::setupMask()
{
	if (mask)
	{
		free(mask);
		mask = NULL;
	}

	mask = (int *) calloc(nn, sizeof(int));
}

void FFT::aboveValueToMask(double value)
{
	setupMask();

	for (int i = 0; i < nn; i++)
	{
		if (data[i][0] >= value)
		{
			mask[i] = 1;
		}
		else
		{
			mask[i] = 0;
		}
	}
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

/*
*    Shift the array in 3D by (nx,ny,nz) pixels
*    Wrap around at the edges
*/
void FFT::shift(long sx, long sy, long sz)
{
	//  printf("Shift: (%li, %li, %li)\n", sx, sy, sz);

	long    x1,y1,z1;
	long    e0,e1;

	fftwf_complex *temp =  (FFTW_DATA_TYPE*) fftwf_malloc(nn*sizeof(FFTW_DATA_TYPE));

	for(long z0=0; z0<nz; z0++)
	{
		z1 = z0 + sz;
		if (z1 < 0) z1 += nz;
		if (z1 >= nz) z1 -= nz;

		for(long y0=0; y0<ny; y0++)
		{
			y1 = y0 + sy;
			if (y1 < 0) y1 += ny;
			if (y1 >= ny) y1 -= ny;

			for(long x0=0; x0<nx; x0++)
			{
				x1 = x0 + sx;
				if(x1 < 0) x1 += nx;
				if(x1 >= nx) x1 -= nx;

				e0 = element(x0,y0,z0);
				e1 = element(x1,y1,z1);


				temp[e0][0] = data[e1][0];
				temp[e0][1] = data[e1][1];
			}
		}
	}

	fftwf_free(data);
	data = temp;
}


void FFT::shiftToCorner(void)
{
	long sx,sy,sz;
	sx = -nx/2;
	sy = -ny/2;
	sz = -nz/2;

	shift(sx,sy,sz);

}

void FFT::shiftToCenter(void)
{
	long sx,sy,sz;
	sx = nx/2;
	sy = ny/2;
	sz = nz/2;

	shift(sx,sy,sz);
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

void FFT::setupBlurring()
{
	_blurAmounts.clear();
	
	for (int i = -1; i < 2; i++)
	{
		for (int j = -1; j < 2; j++)
		{
			for (int k = -1; k < 2; k++)
			{
				int moves = abs(i) + abs(j) + abs(k);
				float factor = normal_distribution(moves, 1.);
				_blurAmounts.push_back(factor);
			}
		}
	}		
	
	_setupBlurring = true;
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
	for (int k = -1; k < 2; k++)
	{
		double sz = z + (double)k;
		for (int j = -1; j < 2; j++)
		{
			double sy = y + (double)j;
			for (int i = -1; i < 2; i++)
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
				
				/*
				lx = (int)lrint(sx);
				ly = (int)lrint(sy);
				lz = (int)lrint(sz);
				
				index = element(lx, ly, lz);
				data[index][0] += factor * real;
				*/

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
							data[index][0] += prop * real * factor;
						}	
					}
				}
			}
		}
	}
}

void FFT::multiplyAll(float value)
{
	for(long i=0; i<nn; i++)
	{
		data[i][0] *= value;
		data[i][1] *= value;
	}
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

	if (fftwf_init_threads() == 0)
	{
		printf("\t\tCould not initialise threads\n");
		exit(1);
	}

	fftwf_plan_with_nthreads(nthreads);

	/*
	*    Read Wisdom from file
	*    Size of fftwf_complex used to determine whether we are using fftwf_ or fftwf_
	*/

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
	_myDims->plan = fftwf_plan_dft_3d((int)nx, (int)ny, (int)nz, data, data, 1, fftw_flags);
	_myDims->iplan = fftwf_plan_dft_3d((int)nx, (int)ny, (int)nz, data, data, -1, fftw_flags);

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
	int tries = 0;

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

double FFT::interpolate(vec3 vox000, size_t im)
{
	vec3 remain = make_vec3(vox000.x - (double)((int)vox000.x),
	                        vox000.y - (double)((int)vox000.y),
	vox000.z - (double)((int)vox000.z));

	long vox000x = vox000.x;
	long vox000y = vox000.y;
	long vox000z = vox000.z;
	long vox000xm = vox000.x + 1;
	long vox000ym = vox000.y + 1;
	long vox000zm = vox000.z + 1;

	collapse(&vox000x, &vox000y, &vox000z);
	collapse(&vox000xm, &vox000ym, &vox000zm);
	
	vox000y  *= nx;
	vox000ym *= nx;
	vox000z  *= nx * ny;
	vox000zm *= nx * ny;

	long int idx000 = vox000x + vox000y + vox000z;
	long int idx100 = vox000xm + vox000y + vox000z;
	long int idx010 = vox000x + vox000ym + vox000z;
	long int idx110 = vox000xm + vox000ym + vox000z;
	long int idx001 = vox000x + vox000y + vox000zm;
	long int idx101 = vox000xm + vox000y + vox000zm;
	long int idx011 = vox000x + vox000ym + vox000zm;
	long int idx111 = vox000xm + vox000ym + vox000zm;

	double val00 = data[idx000][im] * (1 - remain.x) +
	data[idx100][im] * remain.x;
	double val01 = data[idx001][im] * (1 - remain.x) +
	data[idx101][im] * remain.x;
	double val10 = data[idx010][im] * (1 - remain.x) +
	data[idx110][im] * remain.x;
	double val11 = data[idx011][im] * (1 - remain.x) +
	data[idx111][im] * remain.x;

	double val0 = val00 * (1 - remain.y) + val10 * remain.y;
	double val1 = val01 * (1 - remain.y) + val11 * remain.y;

	double value = val0 * (1 - remain.z) + val1 * remain.z;

	return value;
}

double FFT::score(FFTPtr fftCrystal, FFTPtr fftThing, vec3 pos,
                  std::vector<CoordVal> *vals, MapScoreType mapScore)
{
	return operation(fftCrystal, fftThing, pos, mapScore, vals);
}


/*  For multiplying point-wise
*  No assumption that interpolation is not needed.
*/
double FFT::operation(FFTPtr fftEdit, FFTPtr fftConst, vec3 add,
                      MapScoreType mapScoreType, std::vector<CoordVal> *vals,
bool sameScale)
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
	vec3 atomOffset; // in crystal coordinates actually.
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

	vec3 shift = make_vec3((double)(-fftAtom->nx) * 0.5,
	                       (double)(-fftAtom->ny) * 0.5,
	(double)(-fftAtom->nz) * 0.5);
	mat3x3_mult_vec(atomVox2Crystal, &shift);

	/* In crystal voxels at the moment */
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
				if (fftAtom->getReal(atomPos.x, atomPos.y, atomPos.z) <= 10e-6)
				{
					continue;
				}

				/* Find the interpolated value which atomPos falls on */
				double atomReal = 0;

				if (atomPos.x < 0) atomPos.x += fftAtom->nx;
				if (atomPos.y < 0) atomPos.y += fftAtom->ny;
				if (atomPos.z < 0) atomPos.z += fftAtom->nz;

				atomReal = fftAtom->interpolate(atomPos, 0);

				/* We add the atom offset so we don't end up with thousands
				* of atoms at the very centre of our map */
				vec3 finalCrystalVox = vec3_add_vec3(crystalPos, cornerCrystal);

				if (mapScoreType == MapScoreAddNoWrap)
				{
					if (finalCrystalVox.x < -fftCrystal->nx / 2 || 
					    finalCrystalVox.y < -fftCrystal->ny / 2 ||
					    finalCrystalVox.z < -fftCrystal->nz / 2)
					{
						continue;
					}					
					
					if (finalCrystalVox.x > fftCrystal->nx / 2 ||
					    finalCrystalVox.y > fftCrystal->ny / 2 ||
					    finalCrystalVox.z > fftCrystal->nz / 2)
					{
						continue;
					}
				}

				if (finalCrystalVox.x < 0) finalCrystalVox.x += fftCrystal->nx;
				if (finalCrystalVox.y < 0) finalCrystalVox.y += fftCrystal->ny;
				if (finalCrystalVox.z < 0) finalCrystalVox.z += fftCrystal->nz;


				/* Get the index of this final crystal voxel. */
				long crystalIndex = fftCrystal->element(finalCrystalVox.x + 0.5,
				                                        finalCrystalVox.y + 0.5,
				                                        finalCrystalVox.z + 0.5);
				
				count++;

				if (mapScoreType == MapScoreTypeCorrel)
				{
					if ((!fftCrystal->_writeToMaskZero &&
					     fftCrystal->getMask(crystalIndex) == 0))
					{
						continue;
					}

					/* We do NOT need to interpolate */
					double realCryst = fftCrystal->getReal(crystalIndex);

					if (vals)
					{
						CoordVal val;
						val.fo = realCryst;
						val.fc = atomReal;
#ifdef COORDVAL_FULL
						long atomEle = fftAtom->element(atomPos);
						val.pos = atomPos;
						val.mask = 0;
						if (fftAtom->mask)
						{
							val.mask = fftAtom->getMask(atomEle);
						}
#endif

						vals->push_back(val);
					}
				}
				else if (mapScoreType == MapScoreTypeCopyToSmaller)
				{
					double realCryst = fftCrystal->interpolate(finalCrystalVox);
					int ele = fftAtom->element(atomPos.x, atomPos.y,
					                           atomPos.z);


					fftAtom->setElement(ele, realCryst, 0);
				}
				else if (mapScoreType == MapScoreTypeNone ||
				         mapScoreType == MapScoreAddNoWrap)
				{
					/* Add the density to the real value of the crystal voxel.*/

					if (fftCrystal->_writeToMaskZero ||
					    fftCrystal->getMask(crystalIndex) != 0)
					{
						fftCrystal->data[crystalIndex][0] += atomReal * volume;
					}
				}
			}
		}
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

void FFT::printSlice(double zVal)
{
	for (int j = 0; j < ny; j++)
	{
		std::cout << "| ";
		for (int i = 0; i < nx; i++)
		{
			std::string symbol = " ";
			double value = getReal(element(i, j, zVal * nz));

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

vec3 FFT::collapseToRealASU(vec3 frac, CSym::CCP4SPG *spaceGroup)
{
	vec3 ijk = frac;

	for (int l = 0; l < spaceGroup->nsymop; l++)
	{
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

void FFT::applySymmetry(CSym::CCP4SPG *spaceGroup, double maxRes)
{
	fftwf_complex *tempData;
	tempData = (fftwf_complex *)fftwf_malloc(nn * sizeof(FFTW_DATA_TYPE));
	memset(tempData, 0, sizeof(FFTW_DATA_TYPE) * nn);

	int count = 0;
	
//	spaceGroup = CSym::ccp4spg_load_by_ccp4_num(155);
	std::cout << "Applying symmetry operations for space group " << spaceGroup->symbol_xHM;
	std::cout << " (" << spaceGroup->spg_num << ")"  << ": " << std::flush;

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

	
	for (int l = 0; l < spaceGroup->nsymop; l++)
	{
		float *rot = &spaceGroup->invsymop[l].rot[0][0];
		count = 0;

		std::cout << l + 1;
		if (l < spaceGroup->nsymop - 1)
		{
			std::cout << ", ";
		}
		std::cout << std::flush;

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
					
					count++;
					long index = element(i, j, k);
					double myAmp = data[index][0];
					double myPhase = data[index][1];

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
	
	std::cout << "... done." << std::endl;

	memcpy(data, tempData, sizeof(FFTW_DATA_TYPE) * nn);
	free(tempData);
}

void FFT::writeReciprocalToFile(std::string filename, double maxResolution,
                                CSym::CCP4SPG *mtzspg, std::vector<double> unitCell,
mat3x3 real2Frac, FFTPtr data)
{
	double nLimit = nx;
	nLimit = nLimit - ((int)nLimit % 2); // make even
	nLimit /= 2;

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
		unitCell[0] = nx * scales[0];
		unitCell[1] = ny * scales[1];
		unitCell[2] = nz * scales[2];
		unitCell[3] = 90;
		unitCell[4] = 90;
		unitCell[5] = 90;

		mtzspg = CSym::ccp4spg_load_by_ccp4_num(1);
		real2Frac = mat3x3_from_unit_cell(&(unitCell[0]));
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
	for (int i = -nLimit; i < nLimit; i++)
	{
		for (int j = -nLimit; j < nLimit; j++)
		{
			for (int k = -nLimit; k < nLimit; k++)
			{
				bool asu = CSym::ccp4spg_is_in_asu(mtzspg, i, j, k);

				if (!asu)
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
					free = getMask(i, j, k);
				}
				else
				{
					std::cout << "";
				}

				double foAmp = sqrt(foInt);

				// i.e. 0 when mask is free flag.

				/* MTZ file stuff */

				fdata[0] = i;
				fdata[1] = j;
				fdata[2] = k;
				fdata[3] = free;
				fdata[4] = foAmp;
				fdata[5] = calcAmp;
				fdata[6] = 2 * foAmp - calcAmp;
				fdata[7] = phase;
				fdata[8] = foAmp - calcAmp;
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
