//
//  LocalCC.cpp
//  vagabond
//
//  Created by Helen Ginn on 21/10/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#include "fftw3d.h"
#include "LocalCC.h"
#include "maths.h"

double LocalCC::localCorrelation(FFTPtr fft1, FFTPtr fft2)
{
	double sumCC = 0;
	double count = 0;

	for (long k = -fft2->nz / 2; k < fft2->nz / 2; k++)
	{
		for (long j = 0; j < fft2->ny; j++)
		{
			for (long i = 0; i < fft2->nx; i++)
			{
				long base_x = i; long base_y = j; long base_z = k;
				fft2->collapse(&base_x, &base_y, &base_z);
				std::vector<double> val1s, val2s;
				double sumInts = 0;

				for (long z = k - 1; z <= k + 1; z++)
				{
					for (long y = j - 1; y <= j + 1; y++)
					{
						for (long x = i - 1; x <= i + 1; x++)
						{
							long new_x = x; long new_y = y; long new_z = z;
							fft2->collapse(&new_x, &new_y, &new_z);

							if (new_x * base_x < 0 || new_y * base_y < 0 ||
								new_z * base_z < 0)
							{
								continue;
							}

							double int1 = fft1->getIntensity(x, y, z);
							double int2 = fft2->getIntensity(x, y, z);

							val1s.push_back(int1);
							val2s.push_back(int2);

							sumInts += int1 + int2;
						}
					}
				}

				double correl = correlation(val1s, val2s);

				if (correl != correl || sumInts != sumInts) continue;

				sumCC += correl * sumInts;
				count += sumInts;
			}
		}
	}

	return sumCC / count;
}
