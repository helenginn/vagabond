/*
 *  fftw.h
 *  phaseit3d
 *
 *  Created by Anton Barty on 26/07/11.
 *  Copyright 2011 Anton Barty. All rights reserved.
 *
 */

#ifndef fftw3d_h
#define fftw3d_h

#include <fftw3.h>
#include "mat3x3.h"
#include "shared_ptrs.h"

#define FFTW_DATA_TYPE fftwf_complex

#define FFTW_INDEX_LOOP(a) for(long p=0; p<(a).nn; p++)
#define FFTW_SPACE_LOOP(a) for(long z=0;z<(a).nz;z++)for(long y=0;y<(a).ny;y++)for(long x=0;x<(a).nx;x++)

//#define space_loop for(long x=0;x<nx;++x)for(long y=0;y<ny;++y)for(long z=0;z<nz;++z)
#define FFTW_INDEX_LOOP_PRIVATE for(long p=0;p<nn;p++)

#define CPLX_SQR(__data, __p) ( __data[__p][0]*__data[__p][0] + __data[__p][1]*__data[__p][1] )
#define CPLX_ABS(__data, __p) sqrt( __data[__p][0]*__data[__p][0] + __data[__p][1]*__data[__p][1] )
#define CPLX_PHASE(__data, __p) atan2( __data[__p][1], __data[__p][0] )

inline void fftwf_product(fftwf_complex comp1, fftwf_complex comp2, float *result)
{
	result[0] = comp1[0] * comp2[0] - comp1[1] * comp2[1];
	result[1] = 2 * comp1[0] * comp2[1];
}

class FFT {
    
public:
    FFT();
	FFT(FFT &other);
	FFT(long);
    ~FFT();
    
    void create(long);
    void create(long, long, long);

	void setupMask();

	void setMask(long i, MaskType value)
	{
		mask[i] = value;
	}

	MaskType getMask(long i)
	{
		return mask[i];
	}

	long element(long x, long y, long z)
	{
		while (x < 0) x += nx;
		while (x >= nx) x -= nx;

		while (y < 0) y += ny;
		while (y >= ny) y -= ny;

		while (z < 0) z += nz;
		while (z >= nz) z -= nz;

		return x + nx*y + (nx*ny)*z;
	}

	long quickElement(long x, long y, long z)
	{
		return x + nx*y + (nx*ny)*z;
	}

	long elementFromUncorrectedFrac(double xfrac, double yfrac, double zfrac);

    void createFFTWplan(int nthreads=1, int verbose=1, unsigned fftw_flags=FFTW_MEASURE);
    void fft(int direction);
    
    void shift(long, long, long);
    void shiftToCorner(void);
    void shiftToCenter(void);

	double getReal(long index);
	double getReal(long x, long y, long z);

	double getImaginary(long index)
	{
		return data[index][1];
	}

	double getIntensity(long x, long y, long z);
	double getPhase(long x, long y, long z);
	void setReal(double xfrac, double yfrac, double zfrac, double real);
	static void collapseFrac(double *xfrac, double *yfrac, double *zfrac);
    void setAll(float);
    void multiplyAll(float);

    void speedTest(int);

	double interpolate(vec3 fractionalVoxel, bool imaginary = false);

	static void add(FFTPtr fftEdit, FFTPtr fftConst, int scale = 1,
					double addX = 0, double addY = 0, double addZ = 0,
					bool sameScale = false, MaskType type = MaskUnchecked)
	{
		operation(fftEdit, fftConst, scale, addX, addY, addZ, false, type);
	}


	static void multiply(FFTPtr fftEdit, FFTPtr fftConst, int scale = 1,
						 double addX = 0, double addY = 0, double addZ = 0,
						 bool sameScale = false)
	{
		operation(fftEdit, fftConst, scale, addX, addY, addZ, true);
	}

	static void operation(FFTPtr fftEdit, FFTPtr fftConst, int scale = 1,
						  double addX = 0, double addY = 0, double addZ = 0,
						  bool sameScale = false, MaskType type = MaskUnchecked);

	static double score(FFTPtr fftCrystal, FFTPtr fftThing, vec3 position);

	void normalise();

	long int equivalentIndexFor(FFT *other, double realX, double realY, double realZ,
								mat3x3 transform,
								double addX = 0, double addY = 0, double addZ = 0,
								bool sameScale = false);
	long int elementFromFrac(double xFrac, double yFrac, double zFrac);

	inline void setElement(long int index, float real, float imag)
	{
		data[index][0] = real;
		data[index][1] = imag;
	}

	double getScale(int dim)
	{
		return scales[dim];
	}

	void setBasis(mat3x3 mat, double sampleScale);

	void setScales(double val)
	{
		scales[0] = val;
		scales[1] = val;
		scales[2] = val;
		_basis = mat3x3_from_unit_cell(val, val, val, 90., 90., 90.);
		_inverse = mat3x3_inverse(_basis);
	}

	mat3x3 getBasis()
	{
		return _basis;
	}

	mat3x3 getBasisInverse()
	{
		return _inverse;
	}

	void printSlice(bool amplitude = false);

public:
    long nx,ny,nz,nn;
    fftwf_complex *data;
	MaskType *mask; // not char due to cpu speed

	double scales[3];
    
private:
    fftwf_plan plan, iplan;
	bool _made_plan;

	/* Transformation from FFT voxel basis vectors into Angstroms */
	mat3x3 _basis;

	/* Transformation from Angstroms into basis vectors for FFT voxels */
	mat3x3 _inverse;
    
};

#endif /* fftw3d_h */
