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
#include <deque>
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

typedef struct
{
	int nx;
	int ny;
	int nz;
	fftwf_plan plan;
	fftwf_plan iplan;
} FourierDimension;

class FFT {
    
public:
    FFT();
	FFT(FFT &other);
	FFT(long);
    ~FFT();
    
    void create(long);
    void create(long, long, long);

	void setupMask();

	double sumAll()
	{
		double reals = 0;
		for (int i = 0; i < nn; i++)
		{
			reals += data[i][0];
		}

		return reals / (double)nn;
	}

	void setMask(long i, MaskType value)
	{
		mask[i] = value;
	}

	MaskType getMask(long x, long y, long z)
	{
		long elem = element(x, y, z);
		return mask[elem];
	}

	MaskType getMask(long i)
	{
		return mask[i];
	}

	void collapse(long *x, long *y, long *z)
	{
		while (*x < 0) *x += nx;
		while (*x >= nx) *x -= nx;

		while (*y < 0) *y += ny;
		while (*y >= ny) *y -= ny;

		while (*z < 0) *z += nz;
		while (*z >= nz) *z -= nz;
	}

	long element(long x, long y, long z)
	{
		collapse(&x, &y, &z);

		return x + nx*y + (nx*ny)*z;
	}

	long quickElement(long x, long y, long z)
	{
		return x + nx*y + (nx*ny)*z;
	}

	long elementFromUncorrectedFrac(double xfrac, double yfrac, double zfrac);

	void createFFTWplan(int nthreads, unsigned fftw_flags = FFTW_MEASURE);
    void fft(int direction);
    
    void shift(long, long, long);
    void shiftToCorner(void);
    void shiftToCenter(void);
	void shiftToCentre();

	double getReal(long index);
	double getReal(long x, long y, long z);

	double getImaginary(long x, long y, long z);

	double getImaginary(long index)
	{
		return data[index][1];
	}

	double getIntensity(long x, long y, long z);
	double getPhase(long x, long y, long z);
	void setReal(double xfrac, double yfrac, double zfrac, double real);
	void addToReal(double xfrac, double yfrac, double zfrac, double real);
	static void collapseFrac(double *xfrac, double *yfrac, double *zfrac);
	static void multiply(FFTPtr fftEdit, FFTPtr fftConst);
    void setAll(float);
    void multiplyAll(float);

	double interpolate(vec3 fractionalVoxel, bool imaginary = false);

	static void add(FFTPtr fftEdit, FFTPtr fftConst,
					vec3 add, bool sameScale = false)
	{
		operation(fftEdit, fftConst, add);
	}

	static double operation(FFTPtr fftEdit, FFTPtr fftConst, vec3 add,
							bool scoreMe = false, std::vector<double> *xs = NULL,
							std::vector<double> *ys = NULL);

	static double score(FFTPtr fftCrystal, FFTPtr fftThing, vec3 position,
						std::vector<double> *xs = NULL,
						std::vector<double> *ys = NULL);

	void normalise();

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

	void invertScale();

	mat3x3 getBasis()
	{
		return _basis;
	}

	mat3x3 getBasisInverse()
	{
		return _inverse;
	}

	void printSlice(bool amplitude = false);
	void applySymmetry(CSym::CCP4SPG *spaceGroup, bool collapse);

    long nx,ny,nz,nn;
    fftwf_complex *data;
	MaskType *mask; // not char due to cpu speed

	double scales[3];
    
private:
	FourierDimension *_myDims;

	/* Transformation from FFT voxel basis vectors into Angstroms */
	mat3x3 _basis;

	/* Transformation from Angstroms into basis vectors for FFT voxels */
	mat3x3 _inverse;

	static std::deque<FourierDimension> _dimensions;
};

#endif /* fftw3d_h */
