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

#ifndef __vagabond__gradiator__
#define __vagabond__gradiator__

#include "shared_ptrs.h"
#include "vec3.h"
#include <map>

/** \class Gradiator
 *  \brief Class to calculate gradients for individual bonds' effects
 *  on the correlation coefficient with electron density.
 *
 *  Gradiator looks after the voxels covered by a polymer. When asking for
 *  a gradient with the electron density, the gradiator will loop round
 *  all voxel positions to gather change in voxel density with respect
 *  to the whack/kick of a given bond.
 *
 *  Each voxel will have a list of single atom positions (AtomPtr and
 *  element of the ensemble as an integer) which need to be asked for their
 *  respective contributions. This can be given a conservative cutoff to
 *  neighbours only reduce the computational resources.
 **/

typedef struct
{
//	AtomPtr atom;
	int res;
	vec3 pos;
	double mag;
	size_t index;
} SingleAtom;

typedef struct
{
	double change;
	std::vector<vec3> pos;
	std::vector<vec3> dir;
} WhackVal;

typedef struct
{
	vec3 pos; /** Position in crystal */
	double obs; /** observed density value */
	double calc; /** calculated density value */
	std::vector<SingleAtom> nearby_atoms; /* List of atoms to be consulted */
} Voxel;

class Gradiator
{
public:
	Gradiator(PolymerPtr polymer);

	void prepareList();
private:
	double deltaVoxel4Whack(WhackPtr whack, Voxel *vox, int dir);
	double deltaVoxelValue(Voxel *vox);
	double dDistanceTodDensity(AtomPtr a, double curr, double delta);
	double deltaDistance(Voxel *vox, size_t n, WhackPtr whack);
	vec3 whackPosition(Voxel *vox, size_t n, WhackPtr whack);

	double syy();
	double sxy();
	double sxx();
	double sum_y();
	double sum_x();
	double correlationCoefficient();

	void deltaSs4Whack(WhackPtr w, double *dsxx, double *dsxy, int dir);
	double deltaSumX4Whack(WhackPtr w, int dir);
	double deltaCC(WhackPtr w, int dir);
	
	std::vector<Voxel> _voxels;
	std::vector<WhackPtr> _ws;
	std::map<WhackPtr, WhackVal> _whacks;
	std::map<AtomPtr, vec3> _atoms;
	
	PolymerPtr _polymer;
	CrystalPtr _crystal;
};


#endif
