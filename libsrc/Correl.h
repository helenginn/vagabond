// Vagabond
// Copyright (C) 2017-2018 Helen Ginn
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

#ifndef __vagabond__correl__
#define __vagabond__correl__

#include "shared_ptrs.h"
#include "Crystal.h"


class Correl
{
public:
	void setCrystalAndData(CrystalPtr crystal, DiffractionPtr data);

	void localCC();
	void setupShells();
private:
	CrystalPtr _crystal;
	FFTPtr _data;
	DiffractionPtr _diff;
	FFTPtr _fft;

	std::vector<ShellInfo> _shells;
};

#endif
