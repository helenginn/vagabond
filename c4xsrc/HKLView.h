// Cluster4x
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

#ifndef __cluster4x__hklview__
#define __cluster4x__hklview__

#include <hcsrc/mat3x3.h>
#include <libsrc/FFT.h>
#include "ClusterPlot.h"

class KeeperGL;

class HKLView : public ClusterPlot
{
Q_OBJECT
public:
	HKLView(VagFFTPtr fft, double scale);

	virtual void populate();
	virtual void recolour();
	
	virtual std::string axisLabel(int i)
	{
		return "";
	}
	
	virtual size_t axisCount()
	{
		return 3;
	}
private:
	void addPoint(vec3 point, double value, double colour);
	VagFFTPtr _fft;
	double _scale;
};

#endif
