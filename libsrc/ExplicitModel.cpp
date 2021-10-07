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

#include "ExplicitModel.h"
#include "Twist.h"
#include <hcsrc/Fibonacci.h>
#include "Options.h"
#include "CSV.h"
#include "Anisotropicator.h"
#include "Atom.h"

ExplicitModel::ExplicitModel() : Model()
{
	_centre = empty_vec3();
	_lowestZ = FLT_MAX;
	_shift.setZero();
	_rotation.setZero();
	_shift.setParent(this);
	_rotation.setParent(this);
	_changedSamples = true;
	_modifySample = -1;
}

vec3 ExplicitModel::meanOfManyPositions(std::vector<BondSample> *positions)
{
	vec3 sum = make_vec3(0, 0, 0);
	double weights = 0;

	for (int i = 0; i < positions->size(); i++)
	{
		vec3 tmp = (*positions)[i].start;
		vec3_mult(&tmp, (*positions)[i].occupancy);
		sum = vec3_add_vec3(sum, tmp);
		weights += (*positions)[i].occupancy;
	}
	
	if (weights <= 0)
	{
		weights = 1;
	}

	vec3_mult(&sum, 1 / weights);

	return sum;
}

void ExplicitModel::sanityCheck()
{
	if (_storedSamples.size() == 0)
	{
		std::cout << "No stored samples!" << std::endl;
	}

	for (int i = 0; i < _storedSamples.size(); i++)
	{
		vec3 start = _storedSamples[i].start;
		vec3 old_start = _storedSamples[i].old_start;
		mat3x3 basis = _storedSamples[i].basis;
		
		if (start.x != start.x) 
		{
			std::cout << "Start position is nan for " << shortDesc() <<
			std::endl;
			exit(0);
			return;
		}
		if (old_start.x != old_start.x)
		{
			std::cout << "Old start position is nan for " << shortDesc() <<
			std::endl;
			return;
		}
		
		for (int j = 0; j < 9; j++)
		{
			if (basis.vals[j] != basis.vals[j])
			{
//				std::cout << "Basis value " << j << " is nan for " <<
//				shortDesc() << std::endl;
				return;
			}
		}
	}
}

inline vec3 get_rotated_vec(vec3 orig, vec3 centre, mat3x3 rot, vec3 shift)
{
	vec3_subtract_from_vec3(&orig, centre);
	mat3x3_mult_vec(rot, &orig);
	vec3_add_to_vec3(&orig, centre);
	vec3_add_to_vec3(&orig, shift);

	return orig;
}

void ExplicitModel::savePositions()
{
	_savedSamples = _storedSamples;
}

const std::vector<BondSample> &ExplicitModel::getFinalPositions()
{
	if (!_recalcFinal && _storedSamples.size())
	{
		return _storedSamples;	
	}
	
	/* Recalculate final positions */
	
	std::vector<BondSample> *positions = getManyPositions();
	_storedSamples = *positions;
	
	/* Apply twists to each bond. */
	for (int i = 0; i < twistCount() && !isAnchor(); i++)
	{
		TwistPtr twist = _twists[i];
		twist->applyToAnchorSamples(_storedSamples);
	}
	
	/* Add shift if twisting */
	if (twistCount() > 0)
	{
		vec3 shift = _shift.getVec3();
		vec3 r = _rotation.getVec3();
		mat3x3 rot = mat3x3_rotate(r.x, r.y, r.z);

		for (int i = 0; i < _storedSamples.size(); i++)
		{
			vec3 ns = get_rotated_vec(_storedSamples[i].start, 
			                          _centre, rot, shift);
			vec3 os = get_rotated_vec(_storedSamples[i].old_start, 
			                          _centre, rot, shift);
			
			_storedSamples[i].start = ns;
			_storedSamples[i].old_start = os;
			
			mat3x3 tmp = mat3x3_mult_mat3x3(rot, _storedSamples[i].basis);
			_storedSamples[i].basis = tmp;
		}
	}

	std::vector<vec3> posOnly;
	posOnly.reserve(_storedSamples.size());

	for (int i = 0; i < _storedSamples.size(); i++)
	{
		posOnly.push_back(_storedSamples[i].start);
	}
	
	if (_useMutex)
	{
		std::lock_guard<std::mutex> lock(guiLock);

		_absolute = meanOfManyPositions(&_storedSamples);

		if (canFish() || _finalPositions.size() == 0)
		{
			_finalPositions = posOnly;
			_finalFish = _absolute;
		}
	}
	else
	{
		_finalPositions = posOnly;
		_absolute = meanOfManyPositions(&_storedSamples);
	}
	
	/* Deset flag to allow caching */
	
	_recalcFinal = false;

	return _storedSamples;
}

mat3x3 ExplicitModel::makeTorsionBasis(vec3 hPos, vec3 maPos,
                              vec3 miPos, vec3 lPos, double *newAngle)
{
	vec3 a2p = vec3_subtract_vec3(hPos, maPos);
	vec3 a2l = vec3_subtract_vec3(lPos, maPos);
	vec3 a2b = vec3_subtract_vec3(miPos, maPos);

	double sql = vec3_sqlength(a2b);
	double dotp = vec3_dot_vec3(a2p, a2b);
	double dotl = vec3_dot_vec3(a2l, a2b);
	double distp = dotp / sql;
	double distl = dotl / sql;

	vec3 a2bp = a2b;
	vec3 a2bl = a2b;
	vec3_mult(&a2bp, distp);
	vec3_mult(&a2bl, distl);
	vec3 heavy_join = vec3_add_vec3(maPos, a2bp);
	vec3 light_join = vec3_add_vec3(maPos, a2bl);

	vec3 reverse_bond = vec3_subtract_vec3(miPos, maPos);
	vec3 xNew = vec3_subtract_vec3(hPos, heavy_join);
	mat3x3 basis = mat3x3_rhbasis(xNew, reverse_bond);

	if (newAngle)
	{
		vec3 light = vec3_subtract_vec3(lPos, light_join);
		double angle = vec3_angle_with_vec3(light, xNew);

		vec3 recross = vec3_cross_vec3(light, xNew);
		double cosine = vec3_cosine_with_vec3(recross, reverse_bond);

		if (cosine > 0)
		{
			angle *= -1;
		}

		*newAngle = angle;
	}

	return basis;
}

void ExplicitModel::writePositionsToFile(std::string filename)
{
	if (!hasExplicitPositions())
	{
		return;
	}

	CSVPtr csv = CSVPtr(new CSV(4, "x", "y", "z", "k"));
	std::vector<BondSample> positions = getFinalPositions();
	
	for (int i = 0; i < positions.size(); i++)
	{
		csv->addEntry(3, positions[i].start.x, 
		              positions[i].start.y,
		              positions[i].start.z,
		              positions[i].kickValue);
	}

	csv->writeToFile(filename);	
}

void ExplicitModel::refreshPositions()
{
	getFinalPositions();
}

double ExplicitModel::getMeanSquareDeviation()
{
	getAnisotropy(true);
	return _isotropicAverage * 8 * M_PI * M_PI;
}

/* For GUI */
std::vector<vec3> ExplicitModel::fishPositions(vec3 *ave)
{
	std::vector<vec3> positions;
	std::lock_guard<std::mutex> lock(guiLock);
	positions = _finalPositions;
	if (ave)
	{
		*ave = _finalFish;
	}

	return positions;
}

void ExplicitModel::getAnisotropy(bool withKabsch)
{
	if (withKabsch)
	{
		std::vector<BondSample> finals = getFinalPositions();
		std::vector<vec3> finalPoints;

		for (int i = 0; i < finals.size(); i++)
		{
			finalPoints.push_back(finals[i].start);
		}

		if (!finalPoints.size()) return;
		
		Anisotropicator tropicator;
		tropicator.setPoints(finalPoints);
		_realSpaceTensor = tropicator.getTensor();
		_anisotropyExtent = tropicator.anisotropyExtent();
		_smallness = tropicator.smallness();
		_longest = tropicator.longestAxis();
		_isotropicAverage = tropicator.isotropicAverage();
	}
	else
	{
		std::vector<BondSample> *positions = getManyPositions();
		std::vector<vec3> points;

		for (int i = 0; i < positions->size(); i++)
		{
			points.push_back((*positions)[i].start);
		}

		Anisotropicator tropicator;
		tropicator.setPoints(points);
		_smallness = tropicator.smallness();
		_anisotropyExtent = tropicator.anisotropyExtent();
	}
}

vec3 ExplicitModel::longestAxis()
{
	getAnisotropy();
	return _longest;
}

double ExplicitModel::anisotropyExtent(bool withKabsch)
{
	getAnisotropy();
	return _anisotropyExtent;
}

void ExplicitModel::clearTwists()
{
	_twists.clear();
	_shift.setZero();
	_rotation.setZero();

	for (int i = 0; i < twistCount(); i++)
	{
		Twist::setTwist(&*_twists[i], 0.);
	}
	
	propagateChange(30, true);
}

void ExplicitModel::addShifts(RefinementStrategyPtr strategy,
                              AtomGroupPtr clearGroup)
{
	_shift.setCleanup(this, ExplicitModel::propagate);
	_rotation.setCleanup(this, ExplicitModel::propagate);
	_shift.addVecToStrategy(strategy, 0.02, 0.002, "shift");
	return;
	_rotation.addVecToStrategy(strategy, deg2rad(0.2), 
	                           deg2rad(.005), "rotation");
	
}

double ExplicitModel::propagate(void *obj)
{
	static_cast<ExplicitModel *>(obj)->propagateChange(30, true);
	return 0;
}

std::vector<vec3> ExplicitModel::makeCloud(double totalPoints,
                                           double b, double fallOff,
                                           std::vector<double> &occs)
{
	double totalSurfaces = 0;
	double factor = pow(totalPoints, 1./3.) * 2;
	int layers = lrint(factor);
	
	if (layers % 2 == 0)
	{
		layers++;
	}
	
	if (totalPoints < 2)
	{
		layers = 1;
	}
	
	std::vector<double> layerSurfaces;

	/* Work out relative ratios of the surfaces on which points
	 * will be generated. */
	for (int i = 1; i <= layers; i++)
	{
		layerSurfaces.push_back(i * i);
		totalSurfaces += i * i;
	}

	double scale = totalPoints / (double)totalSurfaces;

	double addTotal = 0;
	Fibonacci fib;
	fib.generateLattice(layers, 1);
	std::vector<vec3> directions = fib.getPoints();
	std::vector<vec3> sphereAngles;
	
	if (totalPoints < 2)
	{
		sphereAngles.push_back(empty_vec3());
		return sphereAngles;
	}

	sphereAngles.reserve(totalPoints);
	vec3 yAxis = make_vec3(0, 1, 0);

	for (int j = 0; j < layers; j++)
	{
		vec3 cross = vec3_cross_vec3(directions[j], yAxis);
		
		mat3x3 mat = mat3x3_closest_rot_mat(yAxis, directions[j], cross);

		double frac = (double)(j + 1) / (double)layers;
		double m = sqrt(b) * frac;
		frac *= fallOff;

		int samples = layerSurfaces[j] * scale + 1;
		double offset = 2. / (double)samples;
		
		fib.generateLattice(samples, m);
		std::vector<vec3> points = fib.getPoints();
		
		for (int i = 0; i < points.size(); i++)
		{
			double add = exp(-frac * frac);
			occs.push_back(add);
			addTotal += add;
			
			mat3x3_mult_vec(mat, &points[i]);

			sphereAngles.push_back(points[i]);
		}
	}
	
	for (int i = 0; i < occs.size(); i++)
	{
		occs[i] /= addTotal;
	}

	return sphereAngles;
}

void ExplicitModel::getAverageBasisPos(mat3x3 *aveBasis, 
                                       vec3 *aveStart, 
                                       std::vector<BondSample> *vals)
{
	if (vals == NULL)
	{
		vals = &_storedSamples;
	}

	*aveBasis = make_mat3x3();
	*aveStart = make_vec3(0, 0, 0);

	/* This loop gets average positions for the basis of the bond
	 * and the average start position */
	for (size_t i = 0; i < vals->size(); i++)
	{
		mat3x3_add_mat3x3(aveBasis, vals->at(i).basis);
		*aveStart = vec3_add_vec3(*aveStart, vals->at(i).start);
	}

	double samples = vals->size();
	vec3_mult(aveStart, 1 / samples);
	mat3x3_vectors_to_unity(aveBasis);
}

