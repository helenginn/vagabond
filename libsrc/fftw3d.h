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

#define    FFTW_DATA_TYPE fftwf_complex

#define FFTW_INDEX_LOOP(a) for(long p=0; p<(a).nn; p++)
#define FFTW_SPACE_LOOP(a) for(long z=0;z<(a).nz;z++)for(long y=0;y<(a).ny;y++)for(long x=0;x<(a).nx;x++)

//#define space_loop for(long x=0;x<nx;++x)for(long y=0;y<ny;++y)for(long z=0;z<nz;++z)
#define FFTW_INDEX_LOOP_PRIVATE for(long p=0;p<nn;p++)

#define CPLX_SQR(__data, __p) ( __data[__p][0]*__data[__p][0] + __data[__p][1]*__data[__p][1] )
#define CPLX_ABS(__data, __p) sqrt( __data[__p][0]*__data[__p][0] + __data[__p][1]*__data[__p][1] )
#define CPLX_PHASE(__data, __p) atan2( __data[__p][1], __data[__p][0] )

typedef void (fftwf_operation) (fftwf_complex, fftwf_complex, fftwf_complex*);

void fftwf_product(fftwf_complex comp1, fftwf_complex comp2, fftwf_complex *result);
void fftwf_add(fftwf_complex comp1, fftwf_complex comp2, fftwf_complex *result);


class cFFTW3d {
    
public:
    cFFTW3d();
    cFFTW3d(long);
    ~cFFTW3d();
    
    void create(long);
    void create(long, long, long);
    long element(long, long, long);
    
    void createFFTWplan(int nthreads=1, int verbose=1, unsigned fftw_flags=FFTW_MEASURE);
    void fft(int direction);
    
    void shift(long, long, long);
    void shiftToCorner(void);
    void shiftToCenter(void);

	double getReal(long index);
	double getReal(long x, long y, long z);

	double getIntensity(long x, long y, long z);
	double getPhase(long x, long y, long z);
	void setReal(double xfrac, double yfrac, double zfrac, double real);
	static void collapseFrac(double *xfrac, double *yfrac, double *zfrac);
    void setAll(float);
    void multiplyAll(float);

    void maxreal(void);
    void maxreal(char *);
    void maxreal2(char *);

    void speedTest(int);

	static void add(FFTPtr fftEdit, FFTPtr fftConst, int scale = 1,
					double addX = 0, double addY = 0, double addZ = 0,
					bool sameScale = false)
	{
		operation(fftEdit, fftConst, fftwf_add, scale, addX, addY, addZ);
	}


	static void multiply(FFTPtr fftEdit, FFTPtr fftConst, int scale = 1,
						 double addX = 0, double addY = 0, double addZ = 0,
						 bool sameScale = false)
	{
		operation(fftEdit, fftConst, fftwf_product, scale, addX, addY, addZ, true);
	}

	static void operation(FFTPtr fftEdit, FFTPtr fftConst,
						  fftwf_operation *op, int scale = 1,
						  double addX = 0, double addY = 0, double addZ = 0,
						  bool sameScale = false);

	long int equivalentIndexFor(cFFTW3d *other, double realX, double realY, double realZ,
								mat3x3 transform,
								double addX = 0, double addY = 0, double addZ = 0,
								bool sameScale = false);
	long int elementFromFrac(double xFrac, double yFrac, double zFrac);

	void setElement(long int index, fftwf_complex value)
	{
		data[index][0] = value[0];
		data[index][1] = value[1];
	}

	double getScale(int dim)
	{
		return scales[dim];
	}

	void setMat(mat3x3 mat, double sampleScale);

	void setScales(double val)
	{
		scales[0] = val;
		scales[1] = val;
		scales[2] = val;
		_basis = mat3x3_from_unit_cell(val, val, val, 90., 90., 90.);
		_inverse = mat3x3_inverse(_basis);
	}

	void setSampling(double sampling)
	{
		_sampling = sampling;
	}

	double getSampling()
	{
		return _sampling;
	}

	mat3x3 getBasis()
	{
		return _basis;
	}

	mat3x3 getBasisInverse()
	{
		return _inverse;
	}

	void printSlice();

public:
    long nx,ny,nz,nn;
	double _sampling; // in 1/A
    fftwf_complex *data;

	double scales[3];
    
private:
    fftwf_plan plan, iplan;
	mat3x3 _basis, _inverse;
    
};

#endif /* fftw3d_h */
