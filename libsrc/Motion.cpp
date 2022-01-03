// Vagabond
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

#include "Motion.h"
#include "SpaceSample.h"
#include <hcsrc/Converter.h>
#include "Anchor.h"
#include "Options.h"
#include "Polymer.h"
#include <hcsrc/Fibonacci.h>
#include "Anisotropicator.h"
#include <hcsrc/RefinementNelderMead.h>
#include <hcsrc/RefinementList.h>
#include "FlexLocal.h"

Motion::Motion()
{
	_scale = 1;
	_trans = RefineMat3x3Ptr(new RefineMat3x3(this, NULL));
	_transScale = 1.0;
	_refined = false;
	_allAtoms = AtomGroupPtr(new AtomGroup());
	_allAtoms->setName("motion_all");
	_allBackbone = AtomGroupPtr(new AtomGroup());
	_allBackbone->setName("motion_bb");
	_rotation = new Quat4Refine();
	_r = empty_vec3();
	_displacement = new Quat4Refine();
	_d = empty_vec3();
	_centre = empty_vec3();
}

void Motion::updateAtoms()
{
	_allAtoms->empty();
	_allBackbone->empty();

	for (int i = 0; i < _molecules.size(); i++)
	{
		if (_molecules[i].expired())
		{
			std::cout << "Losing molecule due to deletion elsewhere." << std::endl;
			_molecules.erase(_molecules.begin() + i);
		}

		PolymerPtr pol = ToPolymerPtr(_molecules[i].lock());
		_allAtoms->addAtomsFrom(pol);
		AtomGroupPtr back = pol->getAllBackbone();
		_allBackbone->addAtomsFrom(back);
	}
	
	if (!_refined)
	{
		_centre = _allAtoms->initialCentroid();
	}
	else
	{
//		_centre = _allAtoms->centroid();
	}
}

void Motion::removeAtom(AtomPtr atom)
{
	_allAtoms->removeAtom(atom);
	_allBackbone->removeAtom(atom);
}

void Motion::addToPolymer(PolymerPtr pol)
{
	_molecules.push_back(pol);
	updateAtoms();

	AnchorPtr anch = pol->getAnchorModel();
	anch->addMotion(shared_from_this());
}

void Motion::removeFromPolymer(PolymerPtr pol)
{
	std::cout << "Attempting to remove polymer " << pol->getChainID() 
	<< " from Motion " << getName() << std::endl;

	for (int i = 0; i < _molecules.size(); i++)
	{
		if (!_molecules[i].expired() && _molecules[i].lock() == pol)
		{
			std::cout << "Found it!" << std::endl;
			_molecules.erase(_molecules.begin() + i);
			i--;
		}
		
		if (i < 0)
		{
			break;
		}
	}

	updateAtoms();
}

void Motion::applyTranslations(std::vector<BondSample> &stored,
                               double scale)
{
	mat3x3 translation = _trans->getMat3x3();
	Anisotropicator tropicator;
	tropicator.setTensor(translation);
	mat3x3 trans = tropicator.basis();
	mat3x3 iso = trans;
	
	for (int i = 0; i < 3; i++)
	{
		vec3 col = mat3x3_axis(iso, i);
		vec3_set_length(&col, 1);
		mat3x3_set_axis(&iso, i, col);
	}

	iso = mat3x3_transpose(iso);
	trans = mat3x3_mult_mat3x3(trans, iso);

	vec3 sum_start = empty_vec3();
	vec3 sum_old = empty_vec3();

	for (int i = 0; i < stored.size(); i++)
	{
		vec3_add_to_vec3(&sum_start, stored[i].start);
		vec3_add_to_vec3(&sum_old, stored[i].old_start);
	}
	
	vec3_mult(&sum_start, 1 / (double)stored.size());
	vec3_mult(&sum_old, 1 / (double)stored.size());
	
	vec3 offset = empty_vec3();
	double count = 0;

	for (int i = 0; i < stored.size(); i++)
	{
		vec3 start = stored[i].start;
		vec3 diff = vec3_subtract_vec3(start, sum_start);
		
		if (stored[i].space != NULL && stored[i].space->hasPoints())
		{
			diff = stored[i].space->point3D(i);
			vec3_mult(&diff, scale);
		}


		mat3x3_mult_vec(trans, &diff);
		vec3_mult(&diff, _transScale);
//		vec3 new_start = vec3_add_vec3(sum_start, diff);
//		vec3 old_start = vec3_add_vec3(sum_old, diff);

		offset += diff;
		count++;

		stored[i].start = start + diff;
		stored[i].old_start += diff;
	}
	
	vec3_mult(&offset, 1 / count);

	for (int i = 0; i < stored.size(); i++)
	{
		stored[i].start -= offset;
		stored[i].old_start -= offset;

	}
}

void Motion::addTranslationParameters(RefinementStrategyPtr strategy)
{
	_trans->addTensorToStrategy(strategy, 0.1, 0.0001, "tr");
}

void Motion::addRigidParameters(RefinementStrategyPtr strategy)
{
	_rotation->addVecToStrategy(strategy, deg2rad(1), deg2rad(0.005), 
	                            "rotation");
	_displacement->addVecToStrategy(strategy, 0.01, 0.001, "displacement");

}

void Motion::rigidRefine()
{
	_centre = _allAtoms->centroid();
	*_stream << "\nRefining rigid body: " << _name << std::endl;
	FlexLocal target;
	NelderMeadPtr neld = NelderMeadPtr(new RefinementNelderMead());
	target.setStream(_stream);
	target.attachToStrategy(neld, _allBackbone);
	target.setAtomGroup(_allAtoms);
	neld->setJobName("rigid");
	addRigidParameters(neld);
	neld->setStream(_stream);
	neld->refine();
}

void Motion::refine(bool reciprocal)
{
	_centre = _allAtoms->centroid();
	*_stream << "\nRefining motion: " << _name << std::endl;
	outputStream();
	
	int maxRot = Options::getMaxRotations();
	bool maxed = false;

	FlexLocal target;
	target.setQuick(true);

	Fibonacci fib;
	fib.generateLattice(31, 0.02);
	std::vector<vec3> points = fib.getPoints();

	for (int j = 0; j < maxRot; j++)
	{
		if (j >= librationCount())
		{
			*_stream << std::endl;
			*_stream << "Introducing rotation #" << j << std::endl;
			
			RefinementListPtr list = RefinementListPtr(new RefinementList());
			list->setJobName("rot_search");
			target.setStream(_stream);
			target.attachToStrategy(list, _allBackbone);
			target.recalculateConstant();
			addLibrationParameters(list, j);
			
			std::vector<double> zero = std::vector<double>(3, 0);
			list->addTestSet(zero);
			
			for (int i = 0; i < points.size(); i++)
			{
				std::vector<double> vals;
				vals.push_back(points[i].x);
				vals.push_back(points[i].y);
				vals.push_back(points[i].z);
				list->addTestSet(vals);
			}

			list->refine();
			_allAtoms->refreshPositions();
			target.recalculateConstant();
			list->outputStream();
			
			if (!list->didChange())
			{
				*_stream << "Nevermind, no demand for a rotation." << std::endl;
				deleteLastScrew();
				maxed = true;
			}
			else
			{
				int num = list->getChosen();
				points.erase(points.begin() + num - 1);
			}
		}
		else
		{
//			maxed = true;
		}

		_allAtoms->refreshPositions();
		
		if (maxRot == librationCount() || maxed)
		{
			break;
		}
	}

	for (int i = 0; i < 1; i++)
	{
		NelderMeadPtr neld = NelderMeadPtr(new RefinementNelderMead());
		target.setStream(_stream);
		target.attachToStrategy(neld, _allBackbone);
		target.recalculateConstant();
		neld->setJobName("translation");
		addTranslationParameters(neld);
		
		if (_refined || true)
		{
			addLibrationParameters(neld, -1);
			if (Options::getScrew())
			{
				addScrewParameters(neld, -1);
			}
			neld->setCycles(120);
		}
		
		neld->refine();
		neld->outputStream();
		
		if (!_refined)
		{
			neld->refine();
		}
		_allAtoms->refreshPositions();
	}

	for (int i = 0; i < 1 && false && !_refined; i++)
	{
		NelderMeadPtr neld = NelderMeadPtr(new RefinementNelderMead());
		neld->setJobName("rots_only");
		target.setStream(_stream);
		target.attachToStrategy(neld, _allBackbone);
		target.recalculateConstant();
		neld->setCycles((i + 1) * 50);

		addLibrationParameters(neld, -1);
		
		neld->refine();
		neld->outputStream();
		_allAtoms->refreshPositions();

		if (Options::getScrew())
		{
			NelderMeadPtr neld = NelderMeadPtr(new RefinementNelderMead());
			neld->setJobName("screws_only");
			target.setStream(_stream);
			target.attachToStrategy(neld, _allBackbone);
			target.recalculateConstant();
			neld->setCycles((i + 1) * 50);
			addScrewParameters(neld, -1);
			neld->refine();
			neld->outputStream();
			_allAtoms->refreshPositions();
		}
	}
	
	target.reportTimings();

	_refined = true;
}

void Motion::prepareMotions(std::vector<SpacePoint> &hyperpoints, double scale)
{
	_sOffsets.clear();
	_sRots.clear();
	mat3x3 translation = _trans->getMat3x3();
	Anisotropicator tropicator;
	tropicator.setTensor(translation);
	mat3x3 trans = tropicator.basis();
	mat3x3 iso = trans;
	
	for (int i = 0; i < 3; i++)
	{
		vec3 col = mat3x3_axis(iso, i);
		vec3_set_length(&col, 1);
		mat3x3_set_axis(&iso, i, col);
	}

	iso = mat3x3_transpose(iso);
	trans = mat3x3_mult_mat3x3(trans, iso);

	vec3 offset = empty_vec3();
	double count = 0;

	for (size_t i = 0; i < hyperpoints.size(); i++)
	{
		vec3 diff = make_vec3(hyperpoints[i][0], hyperpoints[i][1],
		                      hyperpoints[i][2]);
		vec3_mult(&diff, scale);
		mat3x3_mult_vec(trans, &diff);
		vec3_mult(&diff, _transScale);
		offset += diff;
		count++;

		_sOffsets.push_back(diff);
	}
	
	offset /= count;

	mat3x3 basis_rot = getOverallRotation();
	mat3x3 all_rots = make_mat3x3();
	
	for (size_t i = 0; i < _sOffsets.size(); i++)
	{
		vec3 diff = _sOffsets[i] + offset;
		mat3x3 rot_only = basis_rot;

		for (int j = 0; j < _quats.size(); j++)
		{
			vec3 quat = _quats[j]->getVec3();
			vec3_mult(&quat, _scale);
			vec3 bVec = make_vec3(1, 0, quat.x / quat.z);
			mat3x3 rotbasis = mat3x3_rhbasis(bVec, quat);
			
			vec3 rot_vec = quat;
			vec3_set_length(&rot_vec, 1);
			double dot = vec3_dot_vec3(diff, quat);
			
			if (rot_vec.x != rot_vec.x || dot != dot)
			{
				continue;
			}

			mat3x3 rot_mat = mat3x3_unit_vec_rotation(rot_vec, dot);
			rot_only = mat3x3_mult_mat3x3(rot_mat, rot_only);
		}

	//	all_rots = mat3x3_mult_mat3x3(rot_only, all_rots);

		diff += _displacement->getVec3();
		_sRots.push_back(rot_only);
		_sOffsets[i] = diff;
	}

	for (size_t i = 0; i < _sRots.size(); i++)
	{
	//	_sRots[i] = mat3x3_mult_mat3x3(reverse, _sRots[i]);
	}
}

void Motion::applyMotions(std::vector<BondSample> &stored, double scale)
{
	if (_sRots.size() == 0)
	{
		vec3 before = empty_vec3();
		for (int i = 0; i < stored.size(); i++)
		{
			vec3_add_to_vec3(&before, stored[i].start);
		}
		vec3_mult(&before, 1 / (double)stored.size());

		applyTranslations(stored, scale);
		applyRotations(stored, scale);

		vec3 after = empty_vec3();
		for (int i = 0; i < stored.size(); i++)
		{
			vec3_add_to_vec3(&after, stored[i].start);
		}
		vec3_mult(&after, 1 / (double)stored.size());

		vec3_subtract_from_vec3(&after, before);

		for (int i = 0; i < stored.size(); i++)
		{
			vec3_add_to_vec3(&stored[i].start, after);
			vec3_add_to_vec3(&stored[i].old_start, after);
		}
		return;
	}
	
	for (size_t i = 0; i < stored.size() && i < _sRots.size(); i++)
	{
		vec3 dir = stored[i].start - _centre;
		mat3x3_mult_vec(_sRots[i], &dir);
		stored[i].basis = mat3x3_mult_mat3x3(_sRots[i], stored[i].basis);
		dir += _centre - _sOffsets[i];
		stored[i].start = dir;
	}
}

void Motion::applyRotations(std::vector<BondSample> &stored, double scale)
{
	vec3 position = empty_vec3();

	for (int i = 0; i < stored.size(); i++)
	{
		vec3_add_to_vec3(&position, stored[i].start);
	}
	
	vec3_mult(&position, 1 / (double)stored.size());
	
	vec3 centre_to_vec = vec3_subtract_vec3(position, _centre);
	mat3x3 basis_rot = getOverallRotation();
	vec3 rot_disp = centre_to_vec;
	mat3x3_mult_vec(basis_rot, &rot_disp);
	vec3_subtract_from_vec3(&rot_disp, centre_to_vec);

	vec3 displacement = _displacement->getVec3();
	vec3_add_to_vec3(&displacement, rot_disp);
	
	vec3 tot_disp = empty_vec3();

	for (int i = 0; i < stored.size(); i++)
	{
		vec3 start = stored[i].start;
		vec3 neg_start = vec3_mult(start, -1);
		vec3 diff = vec3_subtract_vec3(start, position);
		mat3x3 rot_only = basis_rot;

		for (int j = 0; j < _quats.size(); j++)
		{
			vec3 quat = _quats[j]->getVec3();
			vec3_mult(&quat, _scale);
			vec3 bVec = make_vec3(1, 0, quat.x / quat.z);
			mat3x3 rotbasis = mat3x3_rhbasis(bVec, quat);
			
			vec3 screw = _screws[j]->getVec3();
			mat3x3_mult_vec(rotbasis, &screw);
			
			vec3 rot_vec = quat;
			vec3_set_length(&rot_vec, 1);
			double dot = vec3_dot_vec3(diff, quat);
			
			if (rot_vec.x != rot_vec.x || dot != dot)
			{
				continue;
			}

			mat3x3 rot_mat = mat3x3_unit_vec_rotation(rot_vec, dot);
			mat3x3 basis = mat3x3_mult_mat3x3(rot_mat, stored[i].basis); 
			rot_only = mat3x3_mult_mat3x3(rot_mat, rot_only);
			
			vec3_add_to_vec3(&screw, centre_to_vec);
			vec3 shift = screw;
			mat3x3_mult_vec(rot_mat, &shift);
			vec3_subtract_from_vec3(&shift, screw);
			
			if (shift.x != shift.x)
			{
				continue;
			}

			/* just dealing with translation shifts due to rotation! */
			vec3_add_to_vec3(&stored[i].old_start, shift);
			vec3_add_to_vec3(&stored[i].start, shift);
		}
		
		/* changing old_start to match rotation */
		vec3 diff_to_old = vec3_subtract_vec3(stored[i].old_start,
		                                      start);
		mat3x3_mult_vec(rot_only, &diff_to_old);
		vec3 new_old = vec3_add_vec3(start, diff_to_old);
		stored[i].old_start = new_old;

		/* rotate bond */
		mat3x3 basis = mat3x3_mult_mat3x3(rot_only, stored[i].basis);
		stored[i].basis = basis;

		/* displacement from whole-molecule translation */
		vec3_add_to_vec3(&stored[i].old_start, displacement);
		vec3_add_to_vec3(&stored[i].start, displacement);
	}
}

void Motion::addLibrationParameters(RefinementStrategyPtr strategy,
                                      int num)
{
	if (num < 0)
	{
		for (int i = 0; i < _quats.size(); i++)
		{
			std::string rot = "rot" + i_to_str(i);
			_quats[i]->addVecToStrategy(strategy, 0.01, 0.0001, rot);
		}
	}
	else
	{
		if (num >= _quats.size())
		{
			_quats.push_back(new Quat4Refine());
			_screws.push_back(new Quat4Refine());
			num = _quats.size() - 1;
		}

		std::string rot = "rot" + i_to_str(num);
		_quats[num]->addVecToStrategy(strategy, 0.01, 0.0001, rot);
	}
}

void Motion::addScrewParameters(RefinementStrategyPtr strategy,
                                int num)
{
	if (num < 0)
	{
		for (int i = 0; i < _quats.size(); i++)
		{
			std::string screw = "offset" + i_to_str(i);
			_screws[i]->addVec2ToStrategy(strategy, 1.0, 0.05, screw);
		}
	}
	else
	{
		if (num >= _quats.size())
		{
			return;
		}

		std::string screw = "offset" + i_to_str(num);
		_screws[num]->addVec2ToStrategy(strategy, 1.0, 0.05, screw);
	}
}

void Motion::addProperties()
{
	addMat3x3Property("translation", _trans->getMat3x3Ptr());
	
	_tmpQuats.clear();
	_tmpScrews.clear();
	
	for (int i = 0; i < _quats.size(); i++)
	{
		vec3 v = _quats[i]->getVec3();
		vec3 s = _screws[i]->getVec3();
		_tmpQuats.push_back(v);
		_tmpScrews.push_back(s);
	}
	
	addDoubleProperty("scale", &_scale);
	addBoolProperty("refined", &_refined);
	addVec3Property("centre", &_centre);
	addVec3ArrayProperty("rots", &_tmpQuats);
	addVec3ArrayProperty("screws", &_tmpScrews);
	
	addStringProperty("name", &_name);
	
	for (int i = 0; i < _molecules.size(); i++)
	{
		addReference("molecule", _molecules[i].lock());
	}

	_r = _rotation->getVec3();
	_d = _displacement->getVec3();
	addVec3Property("rotation", &_r);
	addVec3Property("displacement", &_d);
}

void Motion::linkReference(BaseParserPtr object, std::string category)
{
	if (category == "molecule")
	{
		MoleculePtr mol = ToMoleculePtr(object);
		_molecules.push_back(mol);
	}
}

void Motion::addObject(ParserPtr object, std::string category)
{
	AtomGroupPtr grp = ToAtomGroupPtr(object);
	if (category == "atoms")
	{
		_allAtoms = grp;
	}
	else if (category == "backbone")
	{
		_allBackbone = grp;
	}
}

void Motion::postParseTidy()
{
	deleteQuats();
	
	for (int i = 0; i < _tmpQuats.size(); i++)
	{
		Quat4Refine *q = new Quat4Refine();
		Quat4Refine *s = new Quat4Refine();
		
		Quat4Refine::setX(q, _tmpQuats[i].x);
		Quat4Refine::setY(q, _tmpQuats[i].y);
		Quat4Refine::setZ(q, _tmpQuats[i].z);

		Quat4Refine::setX(s, _tmpScrews[i].x);
		Quat4Refine::setY(s, _tmpScrews[i].y);
		
		_quats.push_back(q);
		_screws.push_back(s);
	}
	
	_rotation->setVec3(_r);
	_displacement->setVec3(_d);
}

void Motion::deleteLastScrew()
{
	if (_screws.size() == 0)
	{
		return;
	}

	_screws.pop_back();
	_quats.pop_back();
}

void Motion::deleteQuats()
{
	for (int i = 0; i < _quats.size(); i++)
	{
		delete _quats[i];
		_quats[i] = NULL;

		delete _screws[i];
		_screws[i] = NULL;
	}
	
	_quats.clear();
	_screws.clear();
}

Motion::~Motion()
{
	deleteQuats();
}

double Motion::getScale(void *object)
{
	return static_cast<Motion *>(object)->_scale;
}

void Motion::setScale(void *object, double scale)
{
	static_cast<Motion *>(object)->_scale = scale;
}

void Motion::absorbScale()
{
	for (int i = 0; i < _quats.size(); i++)
	{
		vec3 v = _quats[i]->getVec3();
		vec3_mult(&v, _scale);
		_quats[i]->setVec3(v);
	}
	
	_scale = 1;

	mat3x3 mat = _trans->getMat3x3();
	mat3x3_mult_scalar(&mat, _transScale);
	_trans->setMat3x3(mat);

	_transScale = 1;
}

void Motion::reset()
{
	_quats.clear();
	_screws.clear();
	
	mat3x3 cleared = make_mat3x3();
	_trans->setMat3x3(cleared);
}

mat3x3 Motion::getOverallRotation()
{
	vec3 rot = _rotation->getVec3();
	mat3x3 mat = mat3x3_rotate(rot.x, rot.y, rot.z);
	return mat;
}

void Motion::outputStream()
{
	std::ostringstream *o = static_cast<std::ostringstream *>(_stream);
	std::cout << o->str();
	o->str("");
}
