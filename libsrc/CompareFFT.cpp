// vagabond
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

#include "CompareFFT.h"
#include "FFT.h"
#include <iostream>

CompareFFT::CompareFFT()
{
	_setupResolutions = false;
	_resCutoff = -1;
}

void CompareFFT::prepare()
{
	if (!_primary || !_secondary)
	{
		std::cout << "No primary or secondary FFTs set" << std::endl;
		throw -1;
	}
	
	CSym::CCP4SPG *spg = _primary->getSpaceGroup();

	VagFFTPtr other = _secondary;
	if (_tertiary)
	{
		other = _tertiary;
	}

	vec3 nLimits = getNLimits(_primary, other);
	_pairs.clear();

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				FFTPair pair;
				memset(&pair, '\0', sizeof(FFTPair));

				int asu = CSym::ccp4spg_is_in_asu(spg, i, j, k);
				if (!asu)
				{
					pair.skip_score = true;
				}

				if (!_primary->withinBounds(i, j, k))
				{
					pair.skip_score = true;
				}
				
				/* grab indices for each compared FFT */
				long pi = _primary->element(i, j, k);
				long si = _secondary->element(i, j, k);
				long ti = (_tertiary ? _tertiary->element(i, j, k) : 0);
				pair.id1 = pi;
				pair.id2 = si;
				pair.id3 = si;
				
				float real = _primary->getReal(pi);
				float imag = _primary->getImag(pi);
				if (real != real || imag != imag)
				{
					pair.skip_score = true;
				}

				pair.data1[0] = real;
				pair.data1[1] = imag;
				
				if (_tertiary)
				{
					real = _tertiary->getReal(ti);
					imag = _tertiary->getImag(ti);
					if (real != real || imag != imag)
					{
						pair.skip_score = true;
					}

					pair.data3[0] = real;
					pair.data3[1] = imag;
				}

				if (_setupResolutions || _resCutoff > 0)
				{
					double d = _primary->resolution(i, j, k);
					pair.resolution = d;
					
					if (_resCutoff > 0 && d < _resCutoff)
					{
						pair.skip_score = true;
					}
				}

				pair.isFree = _primary->getScratchComponent(pi, 0, 0) < 0.5;
				
				if (pair.isFree)
				{
					pair.skip_score = true;
				}

				_all.push_back(pair);
				
				if (!pair.skip_score)
				{
					_pairs.push_back(pair);
				}
			}
		}
	}
}
