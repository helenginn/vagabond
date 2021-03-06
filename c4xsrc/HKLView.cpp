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

#include "HKLView.h"
#include <libsrc/Crystal.h>
#include <libsrc/FFT.h>
#include "shaders/HKL_fsh.h"
#include "shaders/HKL_vsh.h"

HKLView::HKLView(VagFFTPtr fft, double scale) : ClusterPlot()
{
	_fft = fft;
	_scale = fft->nanlessAverage() * 5;
	
	_fString = hklFsh();
	_vString = hklVsh();
}

void HKLView::addPoint(vec3 point, double value, double colour)
{
	_indices.push_back(_vertices.size());

	Helen3D::Vertex v;
	memset(v.pos, 0, sizeof(Helen3D::Vertex));

	double slide = value / _scale;
	slide = std::min(1., slide);
	v.color[0] = colour / 2;
	v.color[3] = std::max(0., slide);
	
	v.pos[0] = point.x;
	v.pos[1] = point.y;
	v.pos[2] = point.z;

	_vertices.push_back(v);
}

void HKLView::recolour()
{
	return;
}

void HKLView::populate()
{
	mat3x3 uc = _fft->toRecip();
	uc = mat3x3_transpose(uc);
	
	vec3 nLimits = getNLimits(_fft, _fft);

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				long ele = _fft->element(i, j, k);
				double real = _fft->getReal(ele);
				double stdev = _fft->getScratchComponent(ele, 0, 0);
				stdev /= real;
				
				if (real != real)
				{
					continue;
				}

				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(uc, &ijk);
				vec3_mult(&ijk, 10.);
				
				addPoint(ijk, real, stdev);
			}
		}
	}
}
