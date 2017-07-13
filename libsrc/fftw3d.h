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


#define    FFTW_DATA_TYPE fftwf_complex

#define FFTW_INDEX_LOOP(a) for(long p=0; p<(a).nn; p++)
#define FFTW_SPACE_LOOP(a) for(long z=0;z<(a).nz;z++)for(long y=0;y<(a).ny;y++)for(long x=0;x<(a).nx;x++)

//#define space_loop for(long x=0;x<nx;++x)for(long y=0;y<ny;++y)for(long z=0;z<nz;++z)
#define FFTW_INDEX_LOOP_PRIVATE for(long p=0;p<nn;p++)

#define CPLX_SQR(__data, __p) ( __data[__p][0]*__data[__p][0] + __data[__p][1]*__data[__p][1] )
#define CPLX_ABS(__data, __p) sqrt( __data[__p][0]*__data[__p][0] + __data[__p][1]*__data[__p][1] )
#define CPLX_PHASE(__data, __p) atan2( __data[__p][1], __data[__p][0] )



class cFFTW3d {
    
public:
    cFFTW3d();
    cFFTW3d(long);
    ~cFFTW3d();
    
    void create(long);
    void create(long, long, long);
    long element(long, long, long);
    
    
    void createFFTWplan(int nthreads=1, unsigned fftw_flags=FFTW_MEASURE, int verbose=1);	
    void fft(int direction);
    
    void shift(long, long, long);
    void shiftToCorner(void);
    void shiftToCenter(void);

	double getReal(long index);
	double getReal(long x, long y, long z);
	double getIntensity(long x, long y, long z);
	double getPhase(long x, long y, long z);
	void setReal(double xfrac, double yfrac, double zfrac, double real);
    void setAll(float);
    void multiplyAll(float);

    void maxreal(void);
    void maxreal(char *);
    void maxreal2(char *);

    void speedTest(int);
    
	double getScale(int dim)
	{
		return scales[dim];
	}

	void setScale(int dim, double val)
	{
		scales[dim] = val;
	}

	void setSampling(double sampling)
	{
		_sampling = sampling;
	}

	double getSampling()
	{
		return _sampling;
	}

public:
    long nx,ny,nz,nn;
	double _sampling; // in 1/A
    fftwf_complex *data;

	double scales[3];
    
private:
    fftwf_plan plan, iplan;
    
};

#endif /* fftw3d_h */
