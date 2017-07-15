/*
 *  fftw.cpp
 *	A simple object wrapper for 3D FFTs
 *
 *  Created by Anton Barty on 26/07/11.
 *  Copyright 2011 Anton Barty. All rights reserved.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "fftw3d.h"
#include "mat3x3.h"
#include <iostream>

cFFTW3d::cFFTW3d()
{
	nx = 0;
	ny = 0;
	nz = 0;
	nn = 0;
	data = NULL;
}


cFFTW3d::cFFTW3d(long n)
{
	create(n);
}


cFFTW3d::~cFFTW3d() 
{	
	fftwf_free(data); data = NULL;
	fftwf_destroy_plan(plan);
	fftwf_destroy_plan(iplan);
	fftwf_cleanup_threads();
}	

void cFFTW3d::create(long n)
{
	create(n,n,n);
}

void cFFTW3d::create(long nnx, long nny, long nnz)
{
	nx = nnx;
	ny = nny;
	nz = nnz;
	nn = nx*ny*nz;
	free(data); data = NULL;

	data = (FFTW_DATA_TYPE*) fftwf_malloc(nn*sizeof(FFTW_DATA_TYPE));
	if(!data)
	{
		printf("ERROR in fftwData: Malloc failed\n");
		exit(1); 
	}

	memset(data, 0, sizeof(FFTW_DATA_TYPE) * nn);
}


/*
 *	Convert 3D (xyz) triple into 1D array index
 *	x is fastest axis, z is slowest axis
 */
long cFFTW3d::element(long x, long y, long z)
{
	while (x < 0) x += nx;
	while (x >= nx) x -= nx;

	while (y < 0) y += ny;
	while (y >= ny) y -= ny;

	while (z < 0) z += nz;
	while (z >= nz) z -= nz;

	return x + nx*y + (nx*ny)*z;
}


/*
 *	Shift the array in 3D by (nx,ny,nz) pixels
 *	Wrap around at the edges
 */
void cFFTW3d::shift(long sx, long sy, long sz)
{    
    printf("Shift: (%li, %li, %li)\n", sx, sy, sz);
	
	long	x1,y1,z1;
	long	e0,e1;
	
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


void cFFTW3d::shiftToCorner(void)
{
    long sx,sy,sz;
    sx = -nx/2;
    sy = -ny/2;
    sz = -nz/2;
    
    shift(sx,sy,sz);
    
}

void cFFTW3d::shiftToCenter(void)
{
    long sx,sy,sz;
    sx = nx/2;
    sy = ny/2;
    sz = nz/2;
    
    shift(sx,sy,sz);
}

void cFFTW3d::setAll(float value)
{
    for(long i=0; i<nn; i++) 
    {
        data[i][0] = value;
        data[i][1] = value;
    }
}

double cFFTW3d::getPhase(long x, long y, long z)
{
	long index = element(x, y, z);

	double degrees = atan2(data[index][1], data[index][0]) * 180 / M_PI;

	while (degrees >= 360) degrees -= 360;

	while (degrees < 0) degrees += 360;

	return degrees;
}

double cFFTW3d::getIntensity(long x, long y, long z)
{
	long index = element(x, y, z);

	return sqrt(data[index][0] * data[index][0] + data[index][1] * data[index][1]);
}

double cFFTW3d::getReal(long x, long y, long z)
{
	long index = element(x, y, z);

	return data[index][0];
}

double cFFTW3d::getReal(long index)
{
	return data[index][0];
}

void cFFTW3d::setReal(double xfrac, double yfrac, double zfrac, double real)
{
	while (xfrac < 0) xfrac += 1;
	while (xfrac >= 1) xfrac -= 1;

	while (yfrac < 0) yfrac += 1;
	while (yfrac >= 1) yfrac -= 1;

	while (zfrac < 0) zfrac += 1;
	while (zfrac >= 1) zfrac -= 1;
	
	double x = xfrac * nx;
	double y = yfrac * ny;
	double z = zfrac * nz;

	long index = element(x + 0.5, y + 0.5, z + 0.5);

	data[index][0] = real;
	data[index][1] = 0;
}

void cFFTW3d::multiplyAll(float value)
{
    for(long i=0; i<nn; i++)
    {
        data[i][0] *= value;
        data[i][1] *= value;
    }
}

void cFFTW3d::createFFTWplan(int nthreads, unsigned fftw_flags, int verbose)
{

	char	wisdomFile[2048];
	FILE	*fp;

	/*	
	 *	Sanity check
	 */
	if (nx<=0 || ny<=0 || nz<=0)
    {
		printf("Illogical FFT dimensions: %li x %li x %li\n", nx,ny,nz);
		exit(1);
	}

	/* 
	 * 	Threads
	 */
    if (verbose)
    {
		printf("\tUsing nthreads = %i\n",nthreads);
		printf("\tFFT dimensions: %li x %li x %li (unit scales %.2f x %.2f x %.2f Å)\n", nx,ny,nz,
			   scales[0], scales[1], scales[2]);
	}

	if (fftwf_init_threads() == 0)
    {
		printf("\t\tCould not initialise threads\n");
		exit(1);
	}

	fftwf_plan_with_nthreads(nthreads);	
	
	/* 
	 *	Read Wisdom from file 
	 *	Size of fftwf_complex used to determine whether we are using fftwf_ or fftwf_ 
	 */
    
	strcpy(wisdomFile,getenv("HOME"));
	if (sizeof(FFTW_DATA_TYPE) == 2*sizeof(float))
    { 
		strcat(wisdomFile,"/.fftw3f_wisdom");
	}
    else if (sizeof(FFTW_DATA_TYPE) == 2*sizeof(double)) 
	{
    	strcat(wisdomFile,"/.fftw3_wisdom");
	}

    if(verbose)
    {
        printf("\tImporting FFTW wisdom from %s\n",wisdomFile);
	}

    fp = fopen(wisdomFile, "r");
	if (fp != NULL)
    {
		if (!fftwf_import_wisdom_from_file(fp) )
        {
			printf("\t\tError reading wisdom!\n");
		}
        fclose(fp); 	/* be sure to close the file! */
	}
	else
    {
		printf("\t\tCould not open FFTW wisdom file %s\n",wisdomFile);
	}
	
	
	/* 
	 *	Generate FFTW plans 
	 */ 
    if(verbose)
    {
        printf("\tCreating forwards plan....\n");
    }

	plan = fftwf_plan_dft_3d((int)nx, (int)ny, (int)nz, data, data, 1, fftw_flags);

    if(verbose)
    {
        printf("\tCreating inverse plan....\n");
	}
    
    iplan = fftwf_plan_dft_3d((int)nx, (int)ny, (int)nz, data, data, -1, fftw_flags);
	
	
	
	/* 
	 *	Export wisdom to file 
	 */
	fp = fopen(wisdomFile, "w");
	if (fp != NULL)
    {
        if(verbose)
        {
            printf("\tExporting accumulated wisdom to file %s\n",wisdomFile);
        }

		fftwf_export_wisdom_to_file(fp);
	}

	fclose(fp);	
}



/*
 *	Do the FFT
 */
void cFFTW3d::fft(int direction)
{
	if(direction == 1)
	{
    	fftwf_execute(plan); 
	}
    else if (direction == -1)
	{
    	fftwf_execute(iplan); 
	}
    else
    {
		printf("Error: Undefined FFTW direction\n");
		exit(1);
	}
}

void cFFTW3d::speedTest(int nit){
    clock_t			start, start2;
    time_t			startt, endt;
    float			dt;
    float           avg_t;
    
    
    printf("\nCalculating timing for %i FFT/iFFT pairs\n",nit);
    start = clock();
    time(&startt);
    
    for(long it=1; it<=nit; it++) {
        FFTW_INDEX_LOOP_PRIVATE {
            data[p][0] = 1.0;
            data[p][1] = 0.0;
        }
        
        start2 = clock();
        fftwf_execute(plan);
        fftwf_execute(iplan);
        
        time(&endt);
        dt = endt-startt;
        //dt = difftime(endt,startt);
        printf("\t%li : %li^3 FFT/iFFT pair took %3.2f sec\n",it,nx, dt/it);
        
        //avg_t = (float)(clock()-start2);
        //printf("\t%i : %i^3 FFT/iFFT pair took %3.2f sec\n",it,nn, (clock()-start2)/(float)CLOCKS_PER_SEC);
    }
    
    time(&endt);
    dt = endt-startt;
    //dt = difftime(endt,startt);
    avg_t = (clock()-start)/ nit;
    avg_t = avg_t / ((float)CLOCKS_PER_SEC);
    printf("%li^3 FFT/iFFT pair average (CPU time): %3.3f sec\n",nx, avg_t);
    printf("%li^3 FFT/iFFT pair average (human time): %3.3f sec\n",nx, dt/(float) nit);
}


/*
 *  Compute and apply phase factor that maximises reality of the current iterate
 *  Allow for passing of a mask to specify region to be integrated
 */
void cFFTW3d::maxreal(void) {
    maxreal(NULL);
}

/*
 *  Maximise reality
 *  (vector sum and rotation maxtrix)
 */
void cFFTW3d::maxreal(char *mask=NULL) {
    
    float		temp_re, temp_im, temp_sq;
    float       sum_re, sum_im;
    float       corr_re, corr_im;
    float		temp1, temp2;
    float       phase_avg, phase_correction;
    
    /*
     *	Average phase factor (weighted by amplitude to avoid averaging poorly defined weak phases)
     *  Avoid 'if mask==NULL' within loop (put this outside the loop)
     *  Compute average phase directly, or sum up the complex vectors (faster).
     *  Let the compiler optimise redundant variables.
     */
    temp1 = 0;
    temp2 = 0;
    sum_re = 0;
    sum_im = 0;
    // Unmasked version (all data)
    if(mask == NULL ) {
        FFTW_INDEX_LOOP_PRIVATE {
            temp_re = data[p][0];
            temp_im = data[p][1];
            temp_sq = temp_re*temp_re + temp_im*temp_im;
            
            // Weighted by amplitude
            //sum_re += temp_re;
            //sum_im += temp_im;
            
            // Weighted by intensity
            sum_re += temp_re*sqrtf(temp_sq);
            sum_im += temp_im*sqrtf(temp_sq);
        }
    }
    // Masked version (average where mask != 0)
    else {
        FFTW_INDEX_LOOP_PRIVATE {
            if(mask[p] != 0) {
                temp_re = data[p][0];
                temp_im = data[p][1];
                temp_sq = temp_re*temp_re + temp_im*temp_im;

                // Weighted by amplitude
                //sum_re += temp_re;
                //sum_im += temp_im;

                // Weighted by intensity
                sum_re += temp_re*sqrtf(temp_sq);
                sum_im += temp_im*sqrtf(temp_sq);
            }
        }
    }
    phase_avg = atan2f(sum_im, sum_re);
    //printf("%f\n",phase_avg);

    
    phase_correction = phase_avg;
    corr_re = cosf(phase_correction);
    corr_im = sinf(phase_correction);
    
    
    /*
     *	Apply phase factor required to obtain maximum reality
     */
    FFTW_INDEX_LOOP_PRIVATE {
        temp_re = data[p][0];
        temp_im = data[p][1];

        // Complex multiplication (no in-loop atan2, sin, cos)
        data[p][0] = temp_re*corr_re + temp_im*corr_im;
        data[p][1] = -temp_re*corr_im + temp_im*corr_re;
    }
}


/*
 *  Maximise reality
 *  (explicit phase calculation)
 */
void cFFTW3d::maxreal2(char *mask=NULL) {
    
    float		temp_re, temp_im;
    float       sum_re, sum_im;
    float       corr_re, corr_im;
    float		temp1, temp2;
    float       phase_avg, phase_correction, phase_new;
    float       phase, intensity, amplitude;
    
    /*
     *	Average phase factor (weighted by amplitude to avoid averaging poorly defined weak phases)
     *  Avoid 'if mask==NULL' within loop (put this outside the loop)
     *  Compute average phase directly
     *  Redundant variables optimised out by the compiler
     */
    temp1 = 0;
    temp2 = 0;
    sum_re = 0;
    sum_im = 0;
    // Unmasked version (all data)
    if(mask == NULL ) {
        FFTW_INDEX_LOOP_PRIVATE {
            temp_re = data[p][0];
            temp_im = data[p][1];
            
            phase = atan2f(temp_im, temp_re);
            intensity = temp_re*temp_re + temp_im*temp_im;
            amplitude = sqrtf(intensity);
            
            temp1 += phase*intensity;
            temp2 += intensity;
        }
    }
    // Masked version (average where mask != 0)
    else {
        FFTW_INDEX_LOOP_PRIVATE {
            if(mask[p] != 0) {
                temp_re = data[p][0];
                temp_im = data[p][1];
                
                phase = atan2f(temp_im, temp_re);
                intensity = temp_re*temp_re + temp_im*temp_im;
                amplitude = sqrtf(intensity);
                temp1 += phase*amplitude;
                temp2 += amplitude;
            }
        }
    }
    phase_avg = temp1 / temp2;
    //printf("%f\n",phase_avg);
    
    
    phase_correction = phase_avg;
    corr_re = cosf(phase_correction);
    corr_im = sinf(phase_correction);
    
    
    /*
     *	Apply phase factor required to obtain maximum reality
     */
    FFTW_INDEX_LOOP_PRIVATE {
        temp_re = data[p][0];
        temp_im = data[p][1];
        
        phase = atan2f(temp_im, temp_re);
        intensity = temp_re*temp_re + temp_im*temp_im;
        amplitude = sqrtf(intensity);
        
        // Use phase angle
        phase_new = phase + phase_correction;
        data[p][0] = amplitude * cosf(phase_new);
        data[p][1] = amplitude * sinf(phase_new);
    }
}

void writePhaToResolution(double res)
{

}