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

#ifndef __vagabond__CompareFFT__
#define __vagabond__CompareFFT__

#include "shared_ptrs.h"

class CompareFFT
{
public:
	typedef struct
	{
		unsigned long id1;
		unsigned long id2;
		unsigned long id3;
		float data1[2];
		float data3[2];
		double resolution;
		bool isFree;
		bool skip_score;

	} FFTPair;

	CompareFFT();
	
	void setPrimaryFFT(VagFFTPtr fft)
	{
		_primary = fft;
	}
	
	void setSecondaryFFT(VagFFTPtr fft)
	{
		_secondary = fft;
	}
	
	void setTertiaryFFT(VagFFTPtr fft)
	{
		_tertiary = fft;
	}
	
	void prepare();
	
	void setResolutionCutoff(double resol)
	{
		_resCutoff = resol;
	}
	
	void setupResolutions(bool calc)
	{
		_setupResolutions = calc;
	}
	
	size_t pairCount()
	{
		return _pairs.size();
	}

	FFTPair &pair(int i)
	{
		return _pairs[i];
	}
	
	size_t allPairCount()
	{
		return _all.size();
	}

	FFTPair &allPair(int i)
	{
		return _all[i];
	}
private:
	VagFFTPtr _primary;
	VagFFTPtr _secondary;
	VagFFTPtr _tertiary;
	
	bool _setupResolutions;
	double _resCutoff;

	std::vector<FFTPair> _pairs;
	std::vector<FFTPair> _all;
};

#endif
