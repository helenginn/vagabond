//
//  Diffraction.cpp
//  vagabond
//
//  Created by Helen Ginn on 22/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Diffraction.h"
#include "Crystal.h"

void Diffraction::copyToFFT(VagFFTPtr vag)
{
	vec3 nLimits = getNLimits(_fft, _fft);
	
	for (int i = 0; i < vag->nn(); i++)
	{
		vag->setComponent(i, 0, nan(" "));
	}

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				long index = _fft->element(i, j, k);
				long vindex = vag->element(i, j, k);

				double amp = _fft->getReal(index);
				amp *= amp;
				vag->setComponent(vindex, 0, amp);
			}
		}
	}
}

