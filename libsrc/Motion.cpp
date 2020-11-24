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
#include "Converter.h"
#include "Anchor.h"
#include "Options.h"
#include "Polymer.h"
#include "Fibonacci.h"
#include "Anisotropicator.h"
#include "RefinementNelderMead.h"
#include "RefinementList.h"
#include "FlexGlobal.h"

Motion::Motion()
{
	_scale = 1;
	_trans = RefineMat3x3Ptr(new RefineMat3x3(this, NULL));
	_transScale = 1;
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
		_centre = _allAtoms->centroid();
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
                               bool isomorphous)
{
	mat3x3 translation = _trans->getMat3x3();
	Anisotropicator tropicator;
	tropicator.setTensor(translation);
	mat3x3 trans = tropicator.basis();
	
	if (isomorphous)
	{
		for (int i = 0; i < 3; i++)
		{
			vec3 col = mat3x3_axis(trans, i);
			vec3_set_length(&col, 1);
			mat3x3_set_axis(&trans, i, col);
		}
	}

	vec3 sum_start = empty_vec3();
	vec3 sum_old = empty_vec3();

	for (int i = 0; i < stored.size(); i++)
	{
		vec3_add_to_vec3(&sum_start, stored[i].start);
		vec3_add_to_vec3(&sum_old, stored[i].old_start);
	}
	
	vec3_mult(&sum_start, 1 / (double)stored.size());
	vec3_mult(&sum_old, 1 / (double)stored.size());

	for (int i = 0; i < stored.size(); i++)
	{
		vec3 start = stored[i].start;
		vec3 diff = vec3_subtract_vec3(start, sum_start);
		mat3x3_mult_vec(trans, &diff);
		vec3_mult(&diff, _transScale);
		vec3 new_start = vec3_add_vec3(sum_start, diff);
		vec3 old_start = vec3_add_vec3(sum_old, diff);

		stored[i].start = new_start;
		stored[i].old_start = old_start;
	}
}

void Motion::addTranslationParameters(RefinementStrategyPtr strategy)
{
	_trans->addTensorToStrategy(strategy, 0.1, 0.0001, "tr");
}

void Motion::attachTargetToRefinement(RefinementStrategyPtr strategy,
                                      FlexGlobal &target,
                                      bool recip)
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	target.setCrystal(crystal);
	target.setReciprocalRefinement(recip);
	
	if (!recip)
	{
		target.setAtomGroup(_allBackbone);
	}
	else
	{
		target.setAtomGroup(_allAtoms);
	}

	strategy->setVerbose(true);
	strategy->setCycles(100);
	target.getWorkspace().filename = "pre_motion";

	strategy->setEvaluationFunction(FlexGlobal::score, &target);
	FlexGlobal::score(&target);
}

void Motion::rigidRefine()
{
	_centre = _allAtoms->centroid();
	std::cout << "\nRefining rigid body: " << _name << std::endl;
	FlexGlobal target;
	NelderMeadPtr neld = NelderMeadPtr(new RefinementNelderMead());
	attachTargetToRefinement(neld, target, false);
	target.setAtomGroup(_allAtoms);
	neld->setJobName("rigid");
	_rotation->addVecToStrategy(neld, deg2rad(4), deg2rad(0.04), "rotation");
	_displacement->addVecToStrategy(neld, 0.1, 0.001, "displacement");
	neld->refine();
}

void Motion::refine(bool reciprocal)
{
	_centre = _allAtoms->centroid();
	std::cout << "\nRefining motion: " << _name << std::endl;
	
	int maxRot = Options::getMaxRotations();
	bool maxed = false;

	FlexGlobal target;

	Fibonacci fib;
	fib.generateLattice(31, 0.02);
	std::vector<vec3> points = fib.getPoints();

	for (int j = 0; j < maxRot; j++)
	{
		if (j >= librationCount())
		{
			std::cout << std::endl;
			std::cout << "Introducing rotation #" << j << std::endl;
			
			RefinementListPtr list = RefinementListPtr(new RefinementList());
			list->setJobName("rot_search");
			attachTargetToRefinement(list, target, reciprocal);
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
			
			if (!list->didChange())
			{
				std::cout << "Nevermind, no demand for a rotation." << std::endl;
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
		attachTargetToRefinement(neld, target, reciprocal);
		target.recalculateConstant();
		neld->setJobName("translation");
		addTranslationParameters(neld);
		
		Converter conv;
		if (false && _refined)
		{
			conv.setCompareFunction(&target, FlexGlobal::compareParams);
			conv.setStrategy(neld);
		}

		neld->refine();
		_allAtoms->refreshPositions();
	}

	for (int i = 0; i < 2; i++)
	{
		NelderMeadPtr neld = NelderMeadPtr(new RefinementNelderMead());
		neld->setJobName("rots_only");
		attachTargetToRefinement(neld, target, reciprocal);
		target.recalculateConstant();
		neld->setCycles((i + 1) * 50);

		addLibrationParameters(neld, -1);
		
		if (i == 1 && Options::getScrew())
		{
			addScrewParameters(neld, -1);
			neld->setJobName("rots_screws");
		}
		
		if (_refined && false)
		{
			neld->setJobName("rots_screws_trans");
			_trans->addTensorToStrategy(neld, 0.05, 0.00001, "tr");
		}
		
		Converter conv;
		if (false && _refined)
		{
			conv.setCompareFunction(&target, FlexGlobal::compareParams);
			conv.setStrategy(neld);
		}
		
		neld->refine();
		_allAtoms->refreshPositions();
	}
	
	target.reportTimings();

	_refined = true;
}

void Motion::applyMotions(std::vector<BondSample> &stored)
{
	vec3 before = empty_vec3();
	for (int i = 0; i < stored.size(); i++)
	{
		vec3_add_to_vec3(&before, stored[i].start);
	}
	vec3_mult(&before, 1 / (double)stored.size());

	applyTranslations(stored);
	applyRotations(stored);

	vec3 after = empty_vec3();
	for (int i = 0; i < stored.size(); i++)
	{
		vec3_add_to_vec3(&after, stored[i].start);
	}
	vec3_mult(&after, 1 / (double)stored.size());
	
	vec3_subtract_from_vec3(&after, before);
	vec3_mult(&after, -1);

	for (int i = 0; i < stored.size(); i++)
	{
		vec3_subtract_from_vec3(&stored[i].start, after);
		vec3_subtract_from_vec3(&stored[i].old_start, after);
	}
}

void Motion::applyRotations(std::vector<BondSample> &stored)
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
