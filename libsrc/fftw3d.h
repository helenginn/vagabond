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
#include "MapScoreWorkspace.h"
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

typedef enum
{
	MapScoreTypeNone,
	MapScoreTypeCorrel,
	MapScoreTypeCorrelCopy,
	MapScoreTypeRadialMagnitude,
	MapScoreTypeCopyToSmaller,
	MapScoreAddNoWrap,
} MapScoreType;

inline void fftwf_product(fftwf_complex comp1, fftwf_complex comp2, float *result)
{
	result[0] = comp1[0] * comp2[0] - comp1[1] * comp2[1];
	result[1] = 2 * comp1[0] * comp2[1];
}


/** \cond SHOW_FOURIER_DIMENSION */

typedef struct
{
	int nx;
	int ny;
	int nz;
	fftwf_plan plan;
	fftwf_plan iplan;
} FourierDimension;

/** \endcond */

/** \cond SHOW_COORD_VAL */

/** \endcond */

/**
 * \class FFT
 * \brief Deals with three-dimensional grids with voxel data (real-space maps,
 * diffraction data etc.) */

class FFT {

public:
	FFT();
	FFT(FFT &other);
	FFT(long);
	~FFT();

	static void cleanupPlans();
	void create(long);
	void create(long, long, long);
	void copyFrom(FFTPtr other);
	void scaleToFFT(FFTPtr other);

	void setupMask();
	
	double sumImag();
	
	double sumAll()
	{
		double reals = 0;
		for (int i = 0; i < nn; i++)
		{
			reals += fabs(data[i][0]);
		}

		return reals;
	}

	/** Get the average of all real values (ignore imaginary component) */
	double averageAll()
	{
		return sumAll() / (double)nn;
	}

	void setMask(long i, MaskType value)
	{
		mask[i] = value;
	}

	int getMask(long x, long y, long z)
	{
		long elem = element(x, y, z);
		return mask[elem];
	}

	int getMask(long i)
	{
		return mask[i];
	}

	/* If not writing to mask = 0, then mask will not apply to the
	* "add" command in FFT::operation. */
	void avoidWriteToMaskZero(bool set = true)
	{
		_writeToMaskZero = !set;
	}

	/** If parameters x, y, z fall outside of the unit cell (0 < n <= 1)
	 * wrap them round until they do */
	void collapse(long *x, long *y, long *z)
	{
		while (*x < 0) *x += nx;
		while (*x >= nx) *x -= nx;

		while (*y < 0) *y += ny;
		while (*y >= ny) *y -= ny;

		while (*z < 0) *z += nz;
		while (*z >= nz) *z -= nz;
	}

	/** Element corresponding to the vec3 xyz
	 * \param xyz vec3 containing number of voxels along each axis */
	long element(vec3 xyz)
	{
		return element(xyz.x, xyz.y, xyz.z);
	}

	/** Element corresponding to the parameters x, y, z = number of
	 * voxels along each axis */
	long element(long x, long y, long z)
	{
		collapse(&x, &y, &z);

		return x + nx*y + (nx*ny)*z;
	}

	/** Element corresponding to the parameters x, y, z = number of
	 * voxels along each axis, assuming that all parameters stay within
	 * the default unit cell dimensions as supplied */
	long quickElement(long x, long y, long z)
	{
		return x + nx*y + (nx*ny)*z;
	}

	/** Element according to the parameters xfrac, yfrac, zfrac where
	 * these describe a fraction of each unit cell dimension. Will collapse
	 * to default unit cell. */
	long elementFromUncorrectedFrac(double xfrac, double yfrac, double zfrac);

	/** Transfer plans from one FFT object to another if already known
	 * that they are compatible (same nx, ny, nz).*/
	void takePlansFrom(FFTPtr other);
	
	/** Creates an FFTW plan and stores it for future use as a static
	 * storage for the class */
	void createFFTWplan(int nthreads, unsigned fftw_flags = FFTW_MEASURE);
	
	/** Fast fourier transform of associated data
	 * \param direction real -> reciprocal = 1; reciprocal -> real = -1 */
	void fft(int direction);

	/** Move (0, 0, 0) to (nx/2, ny/2, nz/2) and apply same transformation
	 * to entire array, wrapping around if necessary. Useful for seeing
	 * an atom without breaking it up at the densest point, for example */
	void shiftToCentre();

	double getReal(long index)
	{
		return data[index][0];
	}

	double getReal(long x, long y, long z)
	{
		long index = element(x, y, z);

		return data[index][0];
	}

	double getImaginary(long x, long y, long z);

	double getImaginary(long index)
	{
		return data[index][1];
	}

	/** Returns the square of the length of the vector described by the real
	 * and imaginary components of a data index. Takes voxel numbers on
	 * each axis. */
	double getIntensity(long x, long y, long z);

	/** Returns the square of the length of the vector described by the real
	 * and imaginary components of a data index. For a given index directly. */
	double getIntensity(long element);
	
	/** Returns the angle with respect to the real axis for a given
	 * real/imaginary data index. Takes voxel numbers on each axis */
	double getPhase(long x, long y, long z);
	
	/** Set the real coordinate of an index derived from fractional values
	 * to a number. No interpolation. */
	void setReal(double xfrac, double yfrac, double zfrac, double real);
	
	void addBlurredToReal(double xfrac, double yfrac, double zfrac, 
	                      double real);
	void blurRealToImaginary(int i, int j, int k, mat3x3 tensor);
	void addInterpolatedToReal(double sx, double sy, double sz, double val);
	void addInterpolatedToFrac(double fx, double fy, double fz, double val);
	void addToReal(double xfrac, double yfrac, double zfrac, double real);
	static void collapseFrac(double *xfrac, double *yfrac, double *zfrac);
	static void multiply(FFTPtr fftEdit, FFTPtr fftConst);
	void setAll(float);
	void cap(float);
	void multiplyAll(float);
	int setTotal(float value);
	void valueMinus(float value);

	double cubic_interpolate(vec3 vox000, size_t im = false);

	static void add(FFTPtr fftEdit, FFTPtr fftConst,
	                vec3 add, bool sameScale)
	{
		operation(fftEdit, fftConst, add, MapScoreTypeNone, NULL, sameScale);
	}

	static double operation(FFTPtr fftEdit, FFTPtr fftConst, vec3 add,
	                        MapScoreType mapScoreType = MapScoreTypeNone,
	                        std::vector<CoordVal> *vals = NULL, 
	                        bool sameScale = false, bool interp = true);

	static void addSimple(FFTPtr fftEdit, FFTPtr fftConst);
	static double score(FFTPtr fftCrystal, FFTPtr fftThing, vec3 position,
	                    std::vector<CoordVal> *vals = NULL,
	MapScoreType mapScore = MapScoreTypeCorrel);

	void normalise();

	int getMaskFromFrac(vec3 frac)
	{
		long ele = elementFromFrac(frac.x, frac.y, frac.z);
		return mask[ele];
	}

	double getRealFromFrac(vec3 frac)
	{
		frac.x *= nx;
		frac.y *= ny;
		frac.z *= nz;
		return cubic_interpolate(frac, 0);
	}

	long int elementFromFrac(double xFrac, double yFrac, double zFrac);
	vec3 fracFromElement(long int element);
	void nxyzFromElement(long int element, long *x, long *y, long *z);

	inline void setElement(long int index, float real, float imag)
	{
		data[index][0] = real;
		data[index][1] = imag;
	}

	double getScale(int dim)
	{
		return scales[dim];
	}

	void setBasis(mat3x3 mat, double sampleScale = 1);

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

	void bittyShrink(double radius, int num);
	void shrink(double radius);
	void findLimitingValues(double xMin, double xMax, double yMin,
	                        double yMax, double zMin, double zMax,
	                        vec3 *minVals, vec3 *maxVals);

	void addToValueAroundPoint(vec3 pos, double radius, double value,
	                           int bitIndex = -1);
	
	void printCinema();
	void printSlice(int zVal = 0);
	void applySymmetry(CSym::CCP4SPG *spaceGroup, bool silent = false);
	static vec3 collapseToRealASU(vec3 frac, CSym::CCP4SPG *spaceGroup,
	                              int *flipped = NULL);

	void writeReciprocalToFile(std::string filename, double maxResolution = 0,
	                           CSym::CCP4SPG *mtzspg = NULL,
	                           std::vector<double> unitCell = 
	                           std::vector<double>(),
	                           mat3x3 real2Frac = make_mat3x3(),
	                           FFTPtr data = FFTPtr(),
	                           std::vector<double> bins = 
	                           std::vector<double>(),
	                           std::vector<double> ampAves = 
	                           std::vector<double>());

	long nx,ny,nz,nn;
	fftwf_complex *data;
	MaskType *mask;

	static vec3 getPositionInAsu(vec3 vec);
	
	void convertMaskToSolvent(int expTotal);

	double scales[3];
private:
	FourierDimension *_myDims;

	/* Transformation from FFT voxel basis vectors into Angstroms */
	mat3x3 _basis;

	/* Transformation from Angstroms into basis vectors for FFT voxels */
	mat3x3 _inverse;

	bool _writeToMaskZero;
	static std::deque<FourierDimension> _dimensions;

	bool _setupBlurring;
	void setupBlurring();
	std::vector<float> _blurAmounts;
};

#endif /* fftw3d_h */
