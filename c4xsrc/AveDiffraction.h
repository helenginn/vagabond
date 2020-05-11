// cluster4x
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

#ifndef __cluster4x__AveDiffraction__
#define __cluster4x__AveDiffraction__

#include <libsrc/shared_ptrs.h>
#include "Average.h"
#include "MtzFFTPtr.h"

class Group;

class AveDiffraction : public Average
{
public:
	AveDiffraction(Group *group, double maxRes);
	virtual ~AveDiffraction();

	virtual void calculate();
	virtual void findIntercorrelations(Group *other, double **svd);
	
	VagFFTPtr getFFT()
	{
		return _fft;
	}
	
	static std::string unitCellDesc(VagFFTPtr fft);
private:
	virtual double findCorrelation(MtzFFTPtr one, MtzFFTPtr two);
	void scaleIndividualMtz(MtzFFTPtr mtz);
	void scaleIndividuals(Group *other);

	double _maxRes;
	Group *_origGroup;

	VagFFTPtr _fft;
};

#endif
