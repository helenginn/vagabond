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

#define FFTW_DATA_TYPE fftwf_complex

#include "../libccp4/csymlib.h"
#include <fftw3.h>
#include "MapScoreWorkspace.h"
#include "shared_ptrs.h"
#include <hcsrc/mat3x3.h>

class FFT;

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
	int nscratch;
	fftwf_plan atom_real_to_recip;
	fftwf_plan atom_recip_to_real;
	fftwf_plan real_to_recip;
	fftwf_plan recip_to_real;
} FFTDim;

class VagFFT : public boost::enable_shared_from_this<VagFFT>
{
public:
	friend class FFT;

	VagFFT(int nx, int ny, int nz, int nele = 0, int scratches = 0);
	VagFFT(VagFFT &fft, int scratch = -1);
	virtual ~VagFFT();
	
	VagFFTPtr subFFT(int x0, int x1, int y0, int y1, int z0, int z1);
	void prepareShortScratch();
	void addElement(ElementPtr ele);
	void copyScratchElementToPosition(ElementPtr ele);
	
	void wipe();
	void makePlans();
	void setStatus(FFTStatus status)
	{
		_status = status;
	}
	void setSpaceGroup(CSym::CCP4SPG *spg)
	{
		_spg = spg;
	}
	
	CSym::CCP4SPG *getSpaceGroup()
	{
		return _spg;
	}
	
	void reindex(mat3x3 reindex, vec3 v = empty_vec3());
	
	void fft(FFTTransform transform);
	void multiplyFinal(float val);
	void multiplyDotty(float val);
	void multiplyAll(float val);

	void prepareAtomSpace();
	void addAtom(AtomPtr atom);
	double resolution(int i, int j, int k);
	
	/** Applies symmetry operations. If topRes > 0 then resolution is
	 * cut and symmetry application only applies to asymmetric unit,
	 * for speed */
	void applySymmetry(bool silent = true, double topRes = -1,
	                   bool average = true);
	void writeToFile(std::string filename, double maxResolution,
	                 VagFFTPtr data = VagFFTPtr(), 
	                 VagFFTPtr diff = VagFFTPtr(), 
	                 VagFFTPtr calc = VagFFTPtr());
	
	double sumReal(int scratch = -1);
	double nanlessAverage();

	double averageAll()
	{
		return sumReal() / (double)nn();
	}
	
	double sigma();
	void meanSigma(double *mean, double *sigma);
	
	double voxelVolume();
	
	double compareReciprocal(VagFFTPtr other, mat3x3 rotation, 
	                         vec3 translation = empty_vec3(), 
	                         double maxRes = -1, bool adjustCentre = false);

	void setAllReal(double val);
	
	/** a, b, c in Angstroms; alpha, beta, gamma in degrees */
	void setUnitCell(std::vector<double> dims);
	
	std::vector<double> getUnitCell()
	{
		return _unitCell;
	}
	
	/* sets the unit cell dimension in Angstroms and updates other
	 * matrices */
	void setScale(double cubeDim);
	
	void drawSlice(int zVal, std::string filename);
	void printSlice(double zVal = -1, double scale = 1);

	static double operation(VagFFTPtr fftCrystal, VagFFTPtr fftAtom,
	                        MapScoreType mapScoreType, 
	                        std::vector<CoordVal> *vals = NULL,
	                        bool sameScale = false,
	                        bool fcOnly = false,
	                        int step = 1);

	/* simple pair-wise addition of final index reals of two FFTs */
	void addSimple(VagFFTPtr v2);
	
	void copyFrom(VagFFTPtr other);
	
	void copyRealToImaginary();

	double compareReciprocalToScratch(int scratch);
	
	void addToOrigin(vec3 orig)
	{
		vec3_add_to_vec3(&_origin, orig);
	}

	void setOrigin(vec3 orig)
	{
		_origin = orig;
	}
	
	vec3 origin()
	{
		return _origin;
	}
	
	void setScratchComponent(int index, int scratch, int comp, double val);

	void addScratchComponent(int index, int scratch, int comp, double val);
	
	double getScratchComponent(int index, int si, int comp)
	{
		long ii = scratchIndex(index, si);
		return _data[ii][comp];
	}
	
	bool withinBounds(int i, int j, int k);

	/** copy contents of final to scratch map (overwrite) */
	void copyToScratch(int scratch);
	/** add contents of final to scratch map (addition) */
	void addToScratch(int scratch);

	/** add contents of scratch back to final map */
	void addScratchBack(int scratch);
	
	void setActiveScratch(int scratch)
	{
		_activeScratch = scratch;
	}

	/** i = index of nn, j = scratch map number */
	long scratchIndex(int i, int j)
	{
		return i * _stride + _nele * 2 + j + 1;
	}
	
	long finalIndex(int i)
	{
		return i * _stride + _nele * 2;
	}
	
	/** i = index of nn, j = element number */
	long dottyIndex(int i, int j)
	{
		return i * _stride + (j * 2);
	}
	
	/** i = index of nn, j = element number */
	long eleIndex(int i, int j)
	{
		return i * _stride + (j * 2) + 1;
	}

	/** Returns the length of the vector described by the real
	 * and imaginary components of a data index. Takes voxel numbers on
	 * each axis. */
	double getAmplitude(long x, long y, long z)
	{
		double ele = element(x, y, z);
		return getAmplitude(ele);
	}

	double getPhase(int x, int y, int z);

	double minValue();
	double getIntensity(long i);
	double getAmplitude(long i);
	double getPhase(long i);

	double sumAmp(int scratch = -1);

	/* retrieves out of final column */
	double getReal(long i)
	{
		long index = finalIndex(i);
		return _data[index][0];
	}

	/* retrieves out of final column */
	void setComponent(long i, bool imag, double val)
	{
		long index = finalIndex(i);
		_data[index][imag] = val;
	}

	/* retrieves out of final column */
	double getComponent(long i, bool imag)
	{
		long index = finalIndex(i);
		return _data[index][imag];
	}

	/* retrieves out of final column */
	void setReal(int i, int j, int k, double val)
	{
		long ii = element(i, j, k);
		long index = finalIndex(ii);
		_data[index][0] = val;
	}

	/* retrieves out of final column */
	void setReal(long i, double val)
	{
		long index = finalIndex(i);
		_data[index][0] = val;
	}

	/* retrieves out of final column */
	void setImag(long i, double val)
	{
		long index = finalIndex(i);
		_data[index][1] = val;
	}

	/* retrieves out of final column */
	double getReal(int i, int j, int k)
	{
		long ii = element(i, j, k);
		long index = finalIndex(ii);
		return _data[index][0];
	}

	/* retrieves out of final column */
	double getImag(int i, int j, int k)
	{
		long ii = element(i, j, k);
		long index = finalIndex(ii);
		return _data[index][1];
	}

	/* retrieves out of final column */
	void setImag(int i, int j, int k, double val)
	{
		long ii = element(i, j, k);
		long index = finalIndex(ii);
		_data[index][1] = val;
	}
	
	/* retrieves out of final column */
	double getImag(long i)
	{
		long index = finalIndex(i);
		return _data[index][1];
	}

	inline void addToReal(long i, float add)
	{
		long index = finalIndex(i);
		_data[index][0] += add;
	}

	inline void addToImag(long i, float add)
	{
		long index = finalIndex(i);
		_data[index][1] += add;
	}
	
	inline void setElement(long int i, float real, float imag)
	{
		long index = finalIndex(i);
		_data[index][0] = real;
		_data[index][1] = imag;
	}
	
	mat3x3 getRealBasis()
	{
		return _realBasis;
	}
	
	mat3x3 getRecipBasis()
	{
		return _recipBasis;
	}
	
	mat3x3 toReal()
	{
		return _toReal;
	}
	
	mat3x3 toRecip()
	{
		return _toRecip;
	}
	
	double getCubicScale()
	{
		return _toReal.vals[0];
	}

	void adjustNs();
	
	int nx()
	{
		return _nx;
	}
	
	int ny()
	{
		return _ny;
	}
	
	int nz()
	{
		return _nz;
	}
	
	int nn()
	{
		return _nn;
	}
	
	double getBFactor()
	{
		return _bFactor;
	}
	
	void setBFactor(double b)
	{
		_bFactor = b;
	}

	long element(long x, long y, long z)
	{
		collapse(&x, &y, &z);

		return x + _nx*y + (_nx*_ny)*z;
	}

	long element(vec3 xyz)
	{
		return element(xyz.x, xyz.y, xyz.z);
	}
	
	void expandToVoxel(vec3 *frac)
	{
		frac->x *= _nx;
		frac->y *= _ny;
		frac->z *= _nz;
	}

	double getCompFromFrac(vec3 frac, int comp)
	{
		frac.x *= _nx;
		frac.y *= _ny;
		frac.z *= _nz;
		return cubic_interpolate(frac, comp);
	}

	vec3 fracFromElement(long int element);
	void printStatus();

	void addScratchElementToFinal();
	void addInterpolatedToReal(double sx, double sy, 
	                           double sz, double val);

	void switchToComplex();
	void switchToAmplitudePhase();
	void flipCentre();
	vec3 fractionalCentroid();
	vec3 densityCentroid();

	void findLimitingValues(double xMin, double xMax, double yMin,
	                        double yMax, double zMin, double zMax,
	                        vec3 *minVals, vec3 *maxVals);

	void shrink(double radius);
	void bittyShrink(double radius, int num);
	vec3 getSymRelatedPosition(vec3 pos, int i);
	void addToValueAroundPoint(vec3 pos, double radius, double value,
	                           int bitIndex = -1);
	void convertMaskToSolvent(int expTotal);

	/** Remove whole numbers and leave remainder between 0 and 1 for each
	 * x, y, z fraction (i.e. remove unit cells) */
	static void collapseFrac(double *xfrac, double *yfrac, double *zfrac)
	{
		while (*xfrac < 0) *xfrac += 1;
		while (*xfrac >= 1) *xfrac -= 1;

		while (*yfrac < 0) *yfrac += 1;
		while (*yfrac >= 1) *yfrac -= 1;

		while (*zfrac < 0) *zfrac += 1;
		while (*zfrac >= 1) *zfrac -= 1;
	}

	virtual double cubic_interpolate(vec3 vox000, size_t im = false);

	static void setThreads(int num)
	{
		_num_threads = num;
	}

	static void quickFFTs()
	{
		_many_flags = FFTW_ESTIMATE;
	}

	static void patientFFTs()
	{
		_many_flags = FFTW_PATIENT;
	}
protected:
	template <typename T>
	void collapse(T *x, T *y, T *z)
	{
		while (*x < 0) *x += _nx;
		while (*x >= _nx) *x -= _nx;

		while (*y < 0) *y += _ny;
		while (*y >= _ny) *y -= _ny;

		while (*z < 0) *z += _nz;
		while (*z >= _nz) *z -= _nz;
	}

	double getAmplitude(ElementPtr ele, int i, int j, int k);
	int whichColumn(ElementPtr ele);
	void calculateAdjustmentVolumes();
	void addExplicitAtom(AtomPtr atom);
	void addImplicitAtom(AtomPtr atom);
	/** returns true if sane */
	bool sanityCheck();

	void addLessSimple(VagFFTPtr v2);

	double populateImplicit(vec3 centre, vec3 maxVals,
	                        mat3x3 tensor, mat3x3 ellipsoid, double scale);

	/** pre-loaded atom distributions converted to real space in final
	 *  column */
	void separateAtomTransform();

	void setupElements(bool wipe = true);

	std::vector<ElementPtr> _elements;
	std::map<ElementPtr, int> _elementMap;
	std::vector<double> _unitCell;
	
	FFTStatus _status;

	double _bFactor;
	int _nx, _ny, _nz;
	long _nn;
	
	/* individual adjustment volumes for the eight voxel segments when
	 * interpolating */
	double _vols[3];
	
	/** number of elements */
	int _nele;
	int _nscratch;
	int _activeScratch;
	int _stride;

	long _total;
	fftwf_complex *_data;
	fftwf_complex *_lastData;
	float *_shortScratch;
	
	static std::vector<FFTDim *> _dimensions;
	FFTDim *_myDims;
	
	/* _realBasis accounting for unit cell size. Large numbers */
	mat3x3 _toReal;
	/* _recipBasis accounting for unit cell size. Small numbers */
	mat3x3 _toRecip;
	mat3x3 _recipTrans;

	/* small numbers; apply to convert voxel dimension to Angstroms */
	mat3x3 _realBasis;
	/* big numbers; apply to convert Angstrom dimensions to voxel */
	mat3x3 _recipBasis;
	
	mat3x3 _voxelToFrac;
	
	/* denoting where the origin is, e.g. 0,0,0 or atom central pos */
	vec3 _origin;
	
	CSym::CCP4SPG *_spg;

	static int _num_threads;
	static unsigned _many_flags;
	bool _lowResMode;
	bool _setMatrices;
};

#endif
