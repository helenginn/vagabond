/*
 *  fftw.cpp
 *	A simple object wrapper for 3D FFTs
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
#include <vector>
#include "maths.h"

inline void fftwf_add(fftwf_complex comp1, fftwf_complex comp2, float *result)
{
	result[0] = comp1[0];
	result[1] = comp1[1];
	result[0] += comp2[0];
	result[1] += comp2[1];
}


FFT::FFT()
{
	nx = 0;
	ny = 0;
	nz = 0;
	nn = 0;
	data = NULL;
	mask = NULL;
	_made_plan = false;
}


FFT::FFT(long n)
{
	create(n);
}


FFT::~FFT() 
{
	if (_made_plan)
	{
		fftwf_destroy_plan(plan);
		fftwf_destroy_plan(iplan);
	}
	
	fftwf_free(data); data = NULL;
	fftwf_cleanup_threads();
}	

void FFT::create(long n)
{
	create(n,n,n);
}

FFT::FFT(FFT &other)
{
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
		mask = (MaskType *)malloc(nn * sizeof(MaskType));
		memcpy(mask, other.mask, nn * sizeof(MaskType));
	}

	_made_plan = false;

	createFFTWplan(1, false);

	_basis = other._basis;
	_inverse = other._inverse;
}


void FFT::create(long nnx, long nny, long nnz)
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

void FFT::setupMask()
{
	mask = (MaskType *) calloc(nn, sizeof(MaskType));
}

/*
 *	Convert 3D (xyz) triple into 1D array index
 *	x is fastest axis, z is slowest axis
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
	collapseFrac(&xfrac, &yfrac, &zfrac);

	double x = xfrac * nx;
	double y = yfrac * ny;
	double z = zfrac * nz;

	long index = element(x + 0.5, y + 0.5, z + 0.5);

	return index;
}

long FFT::elementFromUncorrectedFrac(double xfrac, double yfrac, double zfrac)
{
	collapseFrac(&xfrac, &yfrac, &zfrac);

	double x = xfrac * nx;
	double y = yfrac * ny;
	double z = zfrac * nz;

	long index = element(x, y, z);

	return index;
}

/*
 *	Shift the array in 3D by (nx,ny,nz) pixels
 *	Wrap around at the edges
 */
void FFT::shift(long sx, long sy, long sz)
{    
  //  printf("Shift: (%li, %li, %li)\n", sx, sy, sz);
	
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

void FFT::setAll(float value)
{
    for(long i=0; i<nn; i++) 
    {
        data[i][0] = value;
        data[i][1] = value;
    }
}

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

double FFT::getReal(long x, long y, long z)
{
	long index = element(x, y, z);

	return data[index][0];
}

double FFT::getImaginary(long x, long y, long z)
{
	long index = element(x, y, z);

	return data[index][1];
}

double FFT::getReal(long index)
{
	return data[index][0];
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

void FFT::multiplyAll(float value)
{
    for(long i=0; i<nn; i++)
    {
        data[i][0] *= value;
        data[i][1] *= value;
    }
}

void FFT::createFFTWplan(int nthreads, int verbose, unsigned fftw_flags)
{
	if (_made_plan)
	{
		return;
	}

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

    if (verbose)
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

		fclose(fp);

		fp = fopen(wisdomFile, "r");

		if (fp != NULL)
		{
			if (!fftwf_import_wisdom_from_file(fp) )
			{
			}

			fclose(fp); 	/* be sure to close the file! */
		}
		else
		{
			printf("\t\tCould not open FFTW wisdom file %s\n",wisdomFile);
		}
	}

	_made_plan = true;
}



/*
 *	Do the FFT
 */
void FFT::fft(int direction)
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

void FFT::speedTest(int nit){
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

double FFT::interpolate(vec3 vox000, bool im)
{
	vec3 remain = make_vec3(fmod(vox000.x, 1), fmod(vox000.y, 1),
							fmod(vox000.z, 1));

	long int idx000 = element(vox000.x, vox000.y, vox000.z);

	long int idx100 = element(vox000.x + 1, vox000.y, vox000.z);
	long int idx010 = element(vox000.x, vox000.y + 1, vox000.z);
	long int idx110 = element(vox000.x + 1, vox000.y + 1, vox000.z);
	long int idx001 = element(vox000.x, vox000.y, vox000.z + 1);
	long int idx101 = element(vox000.x + 1, vox000.y, vox000.z + 1);
	long int idx011 = element(vox000.x, vox000.y + 1, vox000.z + 1);
	long int idx111 = element(vox000.x + 1, vox000.y + 1, vox000.z + 1);

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
				  std::vector<double> *xs, std::vector<double> *ys)
{
	return operation(fftCrystal, fftThing, pos, true, xs, ys);
}


/*  For multiplying point-wise
 *
 */
double FFT::operation(FFTPtr fftEdit, FFTPtr fftConst, vec3 add,
					  bool scoreMe, std::vector<double> *xs,
					  std::vector<double> *ys)
{
	/* I rarely comment something so heavily but I will get confused if
	 * I don't, this time, as I can't soak the protocol into the variable
	 * names. Bear in mind the three coordinate systems:
	 * (a) Angstroms
	 * (b) Crystal voxels
	 * (c) Atom voxels */

	FFT *fftAtom = &*fftConst;
	FFT *fftCrystal = &*fftEdit;
	double volume = fftAtom->getScale(0) * fftAtom->getScale(1)
	* fftAtom->getScale(2);

	/* Bring the fractional coordinate of the atom into range 0 < frac <= 1 */
	FFT::collapseFrac(&add.x, &add.y, &add.z);

	/* Multiply by the relative dimensions of the crystal */
	double multX = add.x * fftCrystal->nx;
	double multY = add.y * fftCrystal->ny;
	double multZ = add.z * fftCrystal->nz;

	/* Get the remainder after subtracting a number of whole voxels */
	vec3 atomOffset;
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

	fftAtom->shiftToCenter();

	vec3 shift = make_vec3((double)(-fftAtom->nx) * 0.5,
						   (double)(-fftAtom->ny) * 0.5,
						   (double)(-fftAtom->nz) * 0.5);
	mat3x3_mult_vec(atomVox2Crystal, &shift);
	vec3 shiftRemainder = make_vec3(fmod(shift.x, 1),
									fmod(shift.y, 1),
									fmod(shift.z, 1));
	vec3 wholeShiftOnly = make_vec3(shift.x - shiftRemainder.x,
									shift.y - shiftRemainder.y,
									shift.z - shiftRemainder.z);

	vec3 cornerCrystal = vec3_add_vec3(wholeShiftOnly, atomWholeCoords);
	atomOffset = vec3_subtract_vec3(atomOffset, shiftRemainder);

	/* We loop around these crystal voxel limits now (ss -> ms -> fs).
	 * We also discard any which happen to go over the limits of our atom voxels
	 * which may occur due to the buffer added above. */

	std::vector<double> crystalVals, thingVals, quickVals;
	crystalVals.reserve(fftAtom->nn); // may be over-estimate, nm.
	thingVals.reserve(fftAtom->nn);

	quickVals.reserve(fftAtom->nn);

	if (scoreMe)
	{
		for (int i = 0; i < fftAtom->nn; i++)
		{
			quickVals.push_back(fftAtom->data[i][0]);
		}
	}

	double meanVal = mean(quickVals);
	int count = 0;

	vec3 atomPos = make_vec3(0, 0, 0);
	for (int k = 0; ; k++)
	{
		for (int j = 0; ; j++)
		{
			for (int i = 0; ; i++)
			{
				/* Position currently in voxel coords - change to atom. */
				vec3 crystalPos = make_vec3(i, j, k);
				atomPos = mat3x3_mult_vec(crystal2AtomVox, crystalPos);
				vec3 pos = atomPos;

				if (atomPos.x > fftAtom->nx)
				{
					break;
				}

				/* Now we must find the relative crystal voxel to write this
				 * density value to, given that the atom is wrapped around
				 * the origin (center). This should work regardless of odd/
				 * even dimension lengths. */

				/* We add the tiny offset which resulted from the atom
				 * falling between two voxels, in atomic voxels */
				vec3 offsetPos = vec3_add_vec3(pos, atomOffset);

				/* Find the interpolated value which offsetPos falls on */
				double atomReal = fftAtom->interpolate(offsetPos, 0);
				double atomImag = 0;

				if (!scoreMe)
				{
					atomImag = fftAtom->interpolate(offsetPos, 1);
				}
				else
				{
					if (atomReal < meanVal / 10)
					{
						continue;
					}

				}

				/* We now convert from atom voxels to crystal voxels */
				mat3x3_mult_vec(atomVox2Crystal, &pos);

				/* We add the atom offset so we don't end up with thousands
				 * of atoms at the very centre of our map */
				vec3 nearlyCrystalVox = vec3_add_vec3(crystalPos, cornerCrystal);
				vec3 finalCrystalVox = (nearlyCrystalVox);

				/* Get the index of this final crystal voxel. */
				long crystalIndex = fftCrystal->element(finalCrystalVox.x,
														finalCrystalVox.y,
														finalCrystalVox.z);

				count++;

				if (scoreMe)
				{
					crystalVals.push_back(fftCrystal->getReal(crystalIndex));
					thingVals.push_back(atomReal);
				}
				else
				{
					/* Add the density to the real value of the crystal voxel.*/
					fftCrystal->data[crystalIndex][0] += atomReal * volume;
					fftCrystal->data[crystalIndex][1] += atomImag * volume;
				}
			}

			if (atomPos.y > fftAtom->ny)
			{
				break;
			}
		}

		if (atomPos.z > fftAtom->nz)
		{
			break;
		}
	}

	//std::cout << "Num: " << count << std::endl;

	if (!scoreMe)
	{
		return 0;
	}

	if (xs && ys)
	{
		*xs = crystalVals;
		*ys = thingVals;
	}

	double correl = correlation(crystalVals, thingVals);

	return correl;
}

/*  For multiplying point-wise
 *
 */
 void FFT::multiply(FFTPtr fftEdit, FFTPtr fftConst)
{
	FFT *fftSmall = &*fftConst;
	FFT *fftBig = &*fftEdit;

	double step = 1;

	for (double k = 0; k < fftSmall->nz; k += step)
	{
		for (double j = 0; j < fftSmall->ny; j += step)
		{
			for (double i = 0; i < fftSmall->nx; i += step)
			{
				long int index = fftSmall->quickElement(i, j, k);

				float real, imag;
				real = fftBig->data[index][0] * fftSmall->data[index][0]
				- fftBig->data[index][1] * fftSmall->data[index][1];
				imag = fftBig->data[index][0] * fftSmall->data[index][1]
				+ fftBig->data[index][1] * fftSmall->data[index][0];

				fftEdit->setElement(index, real, imag);
			}
		}
	}
}


/* To find the same equivalent index bearing in mind the change-of-basis.
 * Assuming that both are centred at the origin. In terms of fractional
 * coordinates but this may need changing when use becomes clear.
 */
long int FFT::equivalentIndexFor(FFT *other, double realX, double realY,
									 double realZ, mat3x3 transform, double addX,
									 double addY, double addZ, bool sameScale)
{
	if (realX > (nx - 1) / 2)
		realX -= (double)nx;
	if (realY > (ny - 1) / 2)
		realY -= (double)ny;
	if (realZ > (nz - 1) / 2)
		realZ -= (double)nz;

	vec3 pos = make_vec3(realX, realY, realZ);

	if (!sameScale)
	{
		/* Get this into Angstrom units */
		mat3x3_mult_vec(transform, &pos);
	}

	pos.x += addX; pos.y += addY; pos.z += addZ;

	long int index = other->element(pos.x, pos.y, pos.z);

	return index;
}

void FFT::printSlice(bool amplitude)
{
	for (int j = 0; j < ny; j++)
	{
		std::cout << "| ";
		for (int i = 0; i < nx; i++)
		{
			std::string symbol = " ";
			double value = getReal(element(i, j, 0));

			if (amplitude)
			{
				value = sqrt(getIntensity(i, j, 0));
			}

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

