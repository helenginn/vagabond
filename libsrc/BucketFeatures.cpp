// Vagabond
// Copyright (C) 2019 Helen Ginn
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

#include "BucketFeatures.h"
#include "Fibonacci.h"
#include "Crystal.h"
#include "fftw3d.h"

void BucketFeatures::addSolvent()
{
	BucketPerStrand::addSolvent();
	_solvCopy = FFTPtr(new FFT(*_solvent));
}

void BucketFeatures::postScaleWork()
{
	for (int i = 0; i < 1; i++)
	{
		enhanceFeatures();
	}
}

void BucketFeatures::enhanceFeatures()
{
	int total = 35;
	Fibonacci fib;
	fib.generateLattice(total, 3.1);
	std::vector<vec3> samples3 = fib.getPoints();

	FFTPtr fft = getCrystal()->getFFT();
	fft->fft(-1);

	mat3x3 real2frac = getCrystal()->getReal2Frac();
	mat3x3 f2r = mat3x3_inverse(real2frac);
	
	for (int k = 0; k < samples3.size(); k++)
	{
		mat3x3_mult_vec(real2frac, &samples3[k]);
	}

	_solvent->setAll(0);
	
	std::cout << "Enhancing solvent features.\r" << std::flush;

	double prog = 0;
	double total_change = 0;
	double count = 0;

	for (int k = 0; k < fft->nz; k++)
	{
		for (int j = 0; j < fft->ny; j++)
		{
			for (int i = 0; i < fft->nx; i++)
			{
				int index = fft->element(i, j, k);
				if (_solvCopy->data[index][0] <= 0)
				{
					continue;
				}

				double fi = (double)i / (double)fft->nx;
				double fj = (double)j / (double)fft->ny;
				int fk = (double)k / (double)fft->nz;
				vec3 pos = make_vec3(fi, fj, fk);

				double far_den = 0;

				for (int l = 0; l < samples3.size(); l++)
				{
					vec3 sample = vec3_add_vec3(pos, samples3[l]);
					double far = fft->getRealFromFrac(sample);

					far_den += far;
				}

				double far_ave = far_den / (double)samples3.size();

				double curr = fft->data[index][0];
				double bit = far_ave - curr;
				
				double change = bit * 0.5;
				double old = _solvCopy->data[index][0];
				
				_solvent->data[index][0] += change * old;
				
				total_change += change;
				count += old;
			}
		}
		
		prog = (double)k / fft->nz;

		printf("\rEnhancing solvent features. |");
		
		for (float p = 0; p < prog; p += 0.02)
		{
			printf("=");
		}
		printf(">");

		for (float p = prog; p < 1; p += 0.02)
		{
			printf(" ");
		}
		printf("| ");

		fflush(stdout);
	}
	
	double corr = total_change / count;

	for (int i = 0; i < _solvent->nn; i++)
	{
		fft->data[i][0] += _solvent->data[i][0] - corr;
	}
	
	std::cout << "Done." << std::endl;
	
	fft->fft(1);
	fft->applySymmetry(getCrystal()->getSpaceGroup());
}
