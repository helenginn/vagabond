// force_down
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

#ifndef __force_down__Wilson__
#define __force_down__Wilson__

#include <libsrc/shared_ptrs.h>
#include <libsrc/Crystal.h>

class Wilson
{
public:
	Wilson(VagFFTPtr vag);

	void forceDown();
	
	void setFilename(std::string fn)
	{
		_filename = fn;
	}
private:
	void findStraightB();
	void splitIntoShells();
	void calculateRatios();
	void correctResolutions();
	double getResidual(size_t maxShell);

	VagFFTPtr _fft;
	std::vector<ShellInfo> _shells;
	
	std::vector<double> _xs;
	std::vector<double> _ys;
	
	std::string _filename;
	double _maxRes;
	size_t _lastShell;
	double _gradient;
	double _intercept;
};

#endif
