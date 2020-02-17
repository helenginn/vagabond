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
	_trans = RefineMat3x3Ptr(new RefineMat3x3(this, NULL));
	_refined = false;
	_allAtoms = AtomGroupPtr(new AtomGroup());
	_allAtoms->setName("motion_all");
	_allBackbone = AtomGroupPtr(new AtomGroup());
	_allBackbone->setName("motion_bb");
	_centre = empty_vec3();
}

void Motion::addToPolymer(PolymerPtr pol)
{
	_allAtoms->addAtomsFrom(pol);
	AtomGroupPtr back = pol->getAllBackbone();
	_allBackbone->addAtomsFrom(back);

	AnchorPtr anch = pol->getAnchorModel();
	anch->addMotion(shared_from_this());
	
	if (!_refined)
	{
		_centre = _allAtoms->initialCentroid();
	}
}

void Motion::translateStartPositions(std::vector<BondSample> &stored)
{
	mat3x3 translation = _trans->getMat3x3();
	Anisotropicator tropicator;
	tropicator.setTensor(translation);
	mat3x3 trans = tropicator.basis();

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
		vec3 new_start = vec3_add_vec3(sum_start, diff);
		vec3 old_start = vec3_add_vec3(sum_old, diff);

		stored[i].start = new_start;
		stored[i].old_start = old_start;
	}
}

void Motion::addTranslationParameters(RefinementStrategyPtr strategy)
{
	_refined = true;
	_trans->addTensorToStrategy(strategy, 0.1, 0.001, "tr");
}

void Motion::attachTargetToRefinement(RefinementStrategyPtr strategy,
                                      FlexGlobal &target)
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	target.setAtomGroup(_allBackbone);
	target.setCrystal(crystal);
	strategy->setVerbose(true);
	strategy->setCycles(100);

	strategy->setEvaluationFunction(FlexGlobal::score, &target);
	FlexGlobal::score(&target);
}

void Motion::refine()
{
	int maxRot = Options::getMaxRotations();
	bool maxed = false;

	FlexGlobal target;
	NelderMeadPtr neld = NelderMeadPtr(new RefinementNelderMead());
	attachTargetToRefinement(neld, target);

	neld->setJobName("translation");
	addTranslationParameters(neld);
	neld->refine();

	_allAtoms->refreshPositions();

	for (int j = 0; j < maxRot; j++)
	{
		if (j >= librationCount())
		{
			std::cout << std::endl;
			std::cout << "Introducing rotation #" << j << std::endl;
			
			RefinementListPtr list = RefinementListPtr(new RefinementList());
			list->setJobName("rot_search");
			attachTargetToRefinement(list, target);
			target.recalculateConstant();
			addLibrationParameters(list, j);

			Fibonacci fib;
			fib.generateLattice(31, 0.02);
			std::vector<vec3> points = fib.getPoints();
			
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
			
			if (!list->didChange())
			{
				std::cout << "Nevermind, no demand for a rotation." << std::endl;
				deleteLastScrew();
				maxed = true;
			}
		}
		else
		{
//			maxed = true;
		}

		_allAtoms->refreshPositions();

		{
			NelderMeadPtr neld = NelderMeadPtr(new RefinementNelderMead());
			neld->setJobName("rots_only");
			attachTargetToRefinement(neld, target);
			target.recalculateConstant();

			addLibrationParameters(neld, -1);
			neld->refine();

		}

		_allAtoms->refreshPositions();
		
		if (maxRot == librationCount() || maxed)
		{
			break;
		}
	}

	{
		NelderMeadPtr neld = NelderMeadPtr(new RefinementNelderMead());
		attachTargetToRefinement(neld, target);
		neld->setJobName("rots_and_offsets");

		addLibrationParameters(neld, -1);
		addScrewParameters(neld, -1);
		neld->refine();
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

	for (int i = 0; i < stored.size(); i++)
	{
		vec3 start = stored[i].start;
		vec3 neg_start = vec3_mult(start, -1);
		vec3 diff = vec3_subtract_vec3(start, position);
		mat3x3 rot_only = make_mat3x3();

		for (int j = 0; j < _quats.size(); j++)
		{
			vec3 quat = _quats[j]->getVec3();
			vec3 bVec = make_vec3(1, 0, quat.x / quat.z);
			mat3x3 rotbasis = mat3x3_rhbasis(bVec, quat);
			mat3x3 transbasis = mat3x3_transpose(rotbasis);

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
			mat3x3 basis = mat3x3_mult_mat3x3(rot_mat, 
			                                  stored[i].basis); 
			rot_only = mat3x3_mult_mat3x3(rot_mat, rot_only);
			
			vec3_subtract_from_vec3(&screw, centre_to_vec);
			vec3 shift = screw;
			mat3x3_mult_vec(rot_mat, &shift);
			vec3_subtract_from_vec3(&shift, screw);
			
			if (shift.x != shift.x)
			{
				continue;
			}


			vec3_add_to_vec3(&stored[i].old_start, shift);
			vec3_add_to_vec3(&stored[i].start, shift);
		}
		
		vec3 diff_to_old = vec3_subtract_vec3(stored[i].old_start,
		                                      start);
		mat3x3_mult_vec(rot_only, &diff_to_old);
		vec3 new_old = vec3_add_vec3(start, diff_to_old);
		stored[i].old_start = new_old;
		mat3x3 basis = mat3x3_mult_mat3x3(rot_only, stored[i].basis);
		stored[i].basis = basis;
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
			_screws[i]->addVec2ToStrategy(strategy, 3.0, 0.01, screw);
		}
	}
	else
	{
		if (num >= _quats.size())
		{
			return;
		}

		std::string screw = "offset" + i_to_str(num);
		_screws[num]->addVec2ToStrategy(strategy, 3.0, 0.01, screw);
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
	
	addBoolProperty("refined", &_refined);
	addVec3Property("centre", &_centre);
	addVec3ArrayProperty("rots", &_tmpQuats);
	addVec3ArrayProperty("screws", &_tmpScrews);
	
	addStringProperty("name", &_name);
	addChild("atoms", _allAtoms);
	addChild("backbone", _allBackbone);
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

void Motion::setScale(double scale)
{
	mat3x3 mat = _trans->getMat3x3();
	mat3x3_scale(&mat, 2, 2, 2);
	_trans->setMat3x3(mat);

}
