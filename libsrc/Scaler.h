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

#ifndef __vagabond__Scaler__
#define __vagabond__Scaler__

#include "shared_ptrs.h"

class Scaler
{
public:
	Scaler(CrystalPtr c, DiffractionPtr d);
	
	void findBestProbeRadius();
	void findProteinSampling();
	
	void setSilent(bool silent)
	{
		_silent = silent;
	}

	void setFull(bool full)
	{
		_full = full;
	}
	
	void setPostAnalysis(bool analysis)
	{
		_analysis = analysis;
	}
	
	void setCycleNumber(int num)
	{
		_cycleNum = num;
	}
	
	void run();
private:
	CrystalPtr _crystal;
	DiffractionPtr _data;
	BucketPtr _bucket;
	VagFFTPtr _fft;

	void fullSolventCalculation();
	void prepareExplicitFFT();
	void absoluteScale();
	void writeSolvent();
	void scaleSolvent();
	void scaleToData();
	void addSolvent();
	void makeMaps();
	void analyse();
	void tidyUp();

	bool _full;
	bool _analysis;
	bool _silent;
	int _cycleNum;
};

#endif
