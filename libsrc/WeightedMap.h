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

#ifndef __vagabond__weighted_map__
#define __vagabond__weighted_map__

#include "shared_ptrs.h"
#include "Crystal.h"

class WeightedMap
{
public:
	WeightedMap();
	
	/** Crystal should be provided after FFT and scaling */
	void setCrystalAndData(CrystalPtr crystal, DiffractionPtr data);

	void createWeightedMaps();
private:
	void calculateFiguresOfMerit();
	void writeCalculatedSlice();
	void writeObservedSlice();
	double stdevForReflection(double fobs, double fcalc, 
	                          double sigfobs, double res);
	double phaseDevForWeight(double weight);

	CrystalPtr _crystal;
	DiffractionPtr _data;
	VagFFTPtr _fft;
	VagFFTPtr _difft;
	std::vector<ShellInfo> _shells;
	double _aveFo;
	double _allWeights;
	double _minZ;
	double _maxZ;

	int shellForResolution(double res);
	double oneMap(FFTPtr scratch, int slice, bool diff);
	void createVagaCoefficients();
	void create2FoFcCoefficients(VagFFTPtr copy = VagFFTPtr(),
	                             bool patt = false);
	void writeFile(VagFFTPtr chosen);
};

#endif
