//
//  Anchor.cpp
//  vagabond
//
//  Created by Helen Ginn 2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "RefineMat3x3.h"
#include "Quat4Refine.h"
#include "Anchor.h"
#include "Motion.h"
#include "Fibonacci.h"
#include "Absolute.h"
#include "Anisotropicator.h"
#include "Options.h"
#include "Polymer.h"
#include "Crystal.h"
#include "Whack.h"
#include "Twist.h"
#include "mat4x4.h"
#include <sstream>

void Anchor::initialise()
{
	_bFactor = 0;
	_alpha = 0;
	_beta = 0;
	_gamma = 0;
	_position = empty_vec3();
	_nDir = empty_vec3();
	_cDir = empty_vec3();
	_trans = RefineMat3x3Ptr(new RefineMat3x3(this, cleanup));
	_disableWhacks = false;
}

Anchor::Anchor(AbsolutePtr absolute) : ExplicitModel()
{
	initialise();
	_bFactor = absolute->getBFactor();
	_position = absolute->getAbsolutePosition();
	_molecule = absolute->getMolecule();
	_atom = absolute->getAtom();
}

void Anchor::deleteQuats()
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

Anchor::~Anchor()
{
	deleteQuats();
}

void Anchor::setNeighbouringAtoms(AtomPtr nPre, AtomPtr nAtom, 
                                  AtomPtr cAtom, AtomPtr cPost)
{
	_nAtom = nAtom;
	_cAtom = cAtom;
	
	vec3 myPos = getAtom()->getInitialPosition();
	vec3 nAtomPos = nAtom->getInitialPosition();
	vec3 cAtomPos = cAtom->getInitialPosition();
	vec3 nPrePos = nPre->getInitialPosition();
	vec3 cPostPos = cPost->getInitialPosition();

	_nDir = vec3_subtract_vec3(nAtomPos, myPos);
	_nDir2 = vec3_subtract_vec3(nPrePos, myPos);
	_cDir = vec3_subtract_vec3(cAtomPos, myPos);
	_cDir2 = vec3_subtract_vec3(cPostPos, myPos);
}

Anchor::Anchor() : ExplicitModel()
{
	initialise();
}

AtomPtr Anchor::getOtherAtom(AtomPtr calling)
{
	if (_nAtom.lock() == calling)
	{
		return _cAtom.lock();
	}
	else if (_cAtom.lock() == calling)
	{
		return _nAtom.lock();
	}
	else
	{
		std::cout << "WATCHOUT" << std::endl;
		return AtomPtr();
	}
}

void Anchor::createLayeredSpherePositions()
{
	_storedSamples.clear();
	
	/* B factor isotropic only atm, get mean square displacement in
	 * each dimension. */
	double meanSqDisp = getBFactor() / (8 * M_PI * M_PI);
	meanSqDisp = sqrt(meanSqDisp);

	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	int totalPoints = crystal->getSampleNum();

	double totalSurfaces = 0;
	double factor = pow(totalPoints, 1./3.) * 2;
	int layers = lrint(factor);
	
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
	
	/* Make Fibonacci lattice for each layer */

	Fibonacci fib;
	_sphereAngles.clear();
	_sphereAngles.reserve(totalPoints);

	for (int j = 0; j < layers; j++)
	{
		double m = meanSqDisp * (double)(j + 1) / (double)layers;

		int samples = layerSurfaces[j] * scale + 1;
		double offset = 2. / (double)samples;
		
		fib.generateLattice(samples, m);
		std::vector<vec3> points = fib.getPoints();
		
		_sphereAngles.insert(_sphereAngles.end(), points.begin(), points.end());
	}
}

void Anchor::createStartPositions(Atom *callAtom)
{
	bool isN = (callAtom == &*(_nAtom.lock()));

	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	int totalPoints = crystal->getSampleNum();
	
	/* Get the rotation matrix for alpha, beta, gamma modifications */
	mat3x3 rot = getAnchorRotation();

	/* Want the direction to be the opposite of the calling bond */
	vec3 direction = isN ? _cDir : _nDir;
	vec3 other = isN ? _cDir2 : _nDir2;
	
	mat3x3_mult_vec(rot, &direction);
	mat3x3_mult_vec(rot, &other);
	
	vec3 empty = empty_vec3();
	double occTotal = 0;

	for (size_t i = 0; i < _sphereAngles.size(); i++)
	{
		vec3 usedPoint = _sphereAngles[i];
		if (totalPoints == 0)
		{
			usedPoint = empty;
		}

		vec3 full = vec3_add_vec3(usedPoint, _position);
		vec3 next = vec3_add_vec3(full, direction);
		vec3 prev = vec3_add_vec3(full, other);
	
		double occ = 1;
		occTotal += occ;
		mat3x3 basis = makeTorsionBasis(prev, next, full, empty);

		BondSample sample;
		sample.basis = basis;
		sample.occupancy = occ;
		sample.torsion = 0;
		sample.old_start = prev; // used instead of atom
		sample.start = full;

		_storedSamples.push_back(sample);
	}

	for (size_t i = 0; i < _storedSamples.size(); i++)
	{
		if (_occupancies.size() == _storedSamples.size())
		{
			_storedSamples[i].occupancy = _occupancies[i];
		}
		else
		{
			_storedSamples[i].occupancy /= occTotal;
		}
	}
}

void Anchor::applyQuaternions()
{
	for (int i = 0; i < _storedSamples.size(); i++)
	{
		vec3 start = _storedSamples[i].start;
		vec3 neg_start = vec3_mult(start, -1);
		vec3 diff = vec3_subtract_vec3(start, _position);
		mat3x3 rot_only = make_mat3x3();

		for (int j = 0; j < _quats.size(); j++)
		{
			vec3 quat = _quats[j]->getVec3();
			mat3x3 rotbasis = mat3x3_ortho_axes(quat);
			mat3x3 transbasis = mat3x3_transpose(rotbasis);

			vec3 screw = _screws[j]->getVec3();
			mat3x3_mult_vec(transbasis, &screw);
			
			vec3 rot_vec = quat;
			vec3_set_length(&rot_vec, 1);
//			double ang = vec3_length(quat);
			double dot = vec3_dot_vec3(diff, quat);
//			dot *= ang;
			
			if (rot_vec.x != rot_vec.x || dot != dot)
			{
				continue;
			}

			mat3x3 rot_mat = mat3x3_unit_vec_rotation(rot_vec, dot);
			mat3x3 basis = mat3x3_mult_mat3x3(rot_mat, 
			                                  _storedSamples[i].basis); 
			rot_only = mat3x3_mult_mat3x3(rot_mat, rot_only);
			
			vec3 shift = screw;
			mat3x3_mult_vec(rot_mat, &shift);
			vec3_subtract_from_vec3(&shift, screw);
			
			if (shift.x != shift.x)
			{
				continue;
			}


			vec3_add_to_vec3(&_storedSamples[i].old_start, shift);
			vec3_add_to_vec3(&_storedSamples[i].start, shift);
		}
		
		vec3 diff_to_old = vec3_subtract_vec3(_storedSamples[i].old_start,
		                                      start);
		mat3x3_mult_vec(rot_only, &diff_to_old);
		vec3 new_old = vec3_add_vec3(start, diff_to_old);
		_storedSamples[i].old_start = new_old;
		mat3x3 basis = mat3x3_mult_mat3x3(rot_only, _storedSamples[i].basis);
		_storedSamples[i].basis = basis;
	}
	
}

void Anchor::atLeastOneMotion()
{
	/* if crystal only has one polymer we do not care */
	CrystalPtr crystal = Options::getActiveCrystal();
	
	if (crystal->polymerCount() <= 1)
	{
		return;
	}
	
	MotionPtr mot = MotionPtr(new Motion());
	PolymerPtr pol = ToPolymerPtr(getMolecule());
	crystal->addMotion(mot);

	mot->setName(getAtom()->shortDesc());
	
	if (pol)
	{
		mot->addToPolymer(pol);
	}
	else
	{
		_motions.push_back(mot);
	}
}

mat3x3 Anchor::getAnchorRotation()
{
	mat3x3 rot = mat3x3_rotate(_alpha, _beta, _gamma);
	
	return rot;
}

void Anchor::applyWholeMotions()
{
	for (int i = 0; i < _motions.size(); i++)
	{
		_motions[i]->translateStartPositions(_storedSamples);
		_motions[i]->applyRotations(_storedSamples);
	}
}

void Anchor::translateStartPositions()
{
	mat3x3 translation = _trans->getMat3x3();
	Anisotropicator tropicator;
	tropicator.setTensor(translation);
	mat3x3 trans = tropicator.basis();

	vec3 sum_start = empty_vec3();
	vec3 sum_old = empty_vec3();

	for (int i = 0; i < _storedSamples.size(); i++)
	{
		vec3_add_to_vec3(&sum_start, _storedSamples[i].start);
		vec3_add_to_vec3(&sum_old, _storedSamples[i].old_start);
	}
	
	vec3_mult(&sum_start, 1 / (double)_storedSamples.size());
	vec3_mult(&sum_old, 1 / (double)_storedSamples.size());

	for (int i = 0; i < _storedSamples.size(); i++)
	{
		vec3 start = _storedSamples[i].start;
		vec3 diff = vec3_subtract_vec3(start, sum_start);
		mat3x3_mult_vec(trans, &diff);
		vec3 new_start = vec3_add_vec3(sum_start, diff);
		vec3 old_start = vec3_add_vec3(sum_old, diff);

		_storedSamples[i].start = new_start;
		_storedSamples[i].old_start = old_start;
	}
}


void Anchor::sanityCheck()
{
	if (vec3_length(_nDir) < 1e-6)
	{
		std::cout << "Anchor N-direction is 0" << std::endl;
	}

	if (vec3_length(_cDir) < 1e-6)
	{
		std::cout << "Anchor C-direction is 0" << std::endl;
	}

	if (vec3_length(_nDir2) < 1e-6)
	{
		std::cout << "Anchor N2-direction is 0" << std::endl;
	}

	if (vec3_length(_cDir2) < 1e-6)
	{
		std::cout << "Anchor C2-direction is 0" << std::endl;
	}
	
	if (_position.x != _position.x)
	{
		std::cout << "Position is nan" << std::endl;
	}
	
	ExplicitModel::sanityCheck();
}

std::string Anchor::shortDesc()
{
	std::ostringstream ss;
	ss << "Anchor_" + getAtom()->shortDesc();
	return ss.str();
}

void Anchor::fixCentroid()
{
	vec3 sum = empty_vec3();

	for (int i = 0; i < _storedSamples.size(); i++)
	{
		vec3 start = _storedSamples[i].start;
		vec3_add_to_vec3(&sum, start);
	}
	
	vec3_mult(&sum, 1 / (double)_storedSamples.size());
	vec3 diff = vec3_subtract_vec3(_position, sum);
	
	for (int i = 0; i < _storedSamples.size(); i++)
	{
		vec3_subtract_from_vec3(&_storedSamples[i].start, diff);
	}
}

void Anchor::recalculateWhacks()
{
	/* Disable all the whacks, then reinstate each whack with a new
	 * base set, then re-enable all whacks. */
	
	_disableWhacks = true;

	for (int i = 0; i < whackCount(); i++)
	{
		_whacks[i]->disable();
	}
	
	for (int i = 0; i < twistCount(); i++)
	{
		_twists[i]->disable();
	}
	
	propagateChange(-1, true);
	
	/* magic matrix should be calculated while all whacks are
	 * off */
	for (int i = 0; i < whackCount(); i++)
	{
		BondPtr child = _whacks[i]->getBond()->downstreamBond(0, 0);
		child->calculateMagicMat();
		child->correctTorsionAngles();
	}

	for (int i = 0; i < whackCount(); i++)
	{
		_whacks[i]->saveSamples();
	}

	for (int i = 0; i < twistCount(); i++)
	{
		_twists[i]->saveSamples();
	}

	_disableWhacks = false;

	for (int i = 0; i < whackCount(); i++)
	{
		_whacks[i]->enable();
	}

	for (int i = 0; i < twistCount(); i++)
	{
		_twists[i]->enable();
	}

	propagateChange(-1, true);
}

std::vector<BondSample> *Anchor::getManyPositions(void *caller)
{
	return getManyPositions(caller, false);	
}

void Anchor::forceRefresh()
{
	getManyPositions(&*_nAtom.lock(), true);
	getManyPositions(&*_cAtom.lock(), true);
}

std::vector<BondSample> *Anchor::getManyPositions(void *caller, bool force)
{
	Atom *callAtom = static_cast<Atom *>(caller);
	
	if (!_samples.count(callAtom))
	{
		SamplePair pair;
		pair.changed = true;
		_samples[callAtom] = pair;
	}
	
	if (!_samples[callAtom].changed)
	{
		return &_samples[callAtom].samples;
	}
	
	createLayeredSpherePositions();
	createStartPositions(callAtom);

	/* Check if number of samples has changed for any reason - if so,
	 * initiate re-caching of initial atom positions for each Whack. */

	if (force)
	{
		recalculateWhacks();
	}
	else
	{
		for (int i = 0; i < _whacks.size() && !_disableWhacks; i++)
		{
			bool refresh = _whacks[i]->needsRefresh(_storedSamples);	

			if (refresh)
			{
				recalculateWhacks();
				break;
			}
		}
	}
	
	if (!_disableWhacks)
	{
		applyWholeMotions();
//		translateStartPositions();
		applyQuaternions();

		fixCentroid();

		/* Apply whacks as normal, if we are not re-caching Whacks. */
		for (int i = 0; i < _whacks.size(); i++)
		{
			WhackPtr whack = _whacks[i];
			whack->applyToAnchorSamples(_storedSamples);
		}
	}

	sanityCheck();
	
	_samples[callAtom].samples = _storedSamples;
	_samples[callAtom].changed = false;
	
	return &_samples[callAtom].samples;
}

void Anchor::addProperties()
{
	addDoubleProperty("bfactor", &_bFactor);
	addDoubleProperty("alpha", &_alpha);
	addDoubleProperty("beta", &_beta);
	addDoubleProperty("gamma", &_gamma);
	addVec3Property("position", &_position);
	addReference("atom", _atom.lock());
	addReference("n_atom", _nAtom.lock());
	addVec3Property("n_dir", &_nDir);
	addVec3Property("c_dir", &_cDir);
	addVec3Property("pre_n", &_nDir2);
	addVec3Property("post_c", &_cDir2);
	addMat3x3Property("translation", _trans->getMat3x3Ptr(), true);
	addReference("c_atom", _cAtom.lock());
	
	for (int i = 0; i < whackCount(); i++)
	{
		addChild("whack", _whacks[i]);
	}
	
	for (size_t i = 0; i < motionCount(); i++)
	{
		addReference("motion", _motions[i]);
	}
	
	_tmpQuats.clear();
	_tmpScrews.clear();
	
	for (int i = 0; i < _quats.size(); i++)
	{
		vec3 v = _quats[i]->getVec3();
		vec3 s = _screws[i]->getVec3();
		_tmpQuats.push_back(v);
		_tmpScrews.push_back(s);
	}
	
	addVec3ArrayProperty("rots", &_tmpQuats, true);
	addVec3ArrayProperty("screws", &_tmpScrews, true);
	
	Model::addProperties();
}

void Anchor::postParseTidy()
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

std::string Anchor::getParserIdentifier()
{
	return "anchor_" + i_to_str(getAtom()->getAtomNum());
}

void Anchor::addObject(ParserPtr object, std::string category)
{
	if (category == "whack")
	{
		WhackPtr whack = ToWhackPtr(object);
		addWhack(whack);
	}
}

void Anchor::linkReference(BaseParserPtr object, std::string category)
{
	AtomPtr atom = ToAtomPtr(object);

	if (category == "atom")
	{
		_atom = atom;
	}
	else if (category == "n_atom")
	{
		_nAtom = atom;
	}
	else if (category == "c_atom")
	{
		_cAtom = atom;
	}
	else if (category == "motion")
	{
		MotionPtr motion = ToMotionPtr(object);
		addMotion(motion);
	}
}

void Anchor::propagateChange(int depth, bool refresh)
{
	for (std::map<Atom *, SamplePair>::iterator it = _samples.begin();
	     it != _samples.end(); it++)
	{
		it->second.changed = true;
	}

	_cAtom.lock()->getModel()->propagateChange(depth, refresh);
	_nAtom.lock()->getModel()->propagateChange(depth, refresh);

	/* Will force recalculation of final positions */
	Model::propagateChange(depth, refresh);
	
	if (refresh)
	{
		getFinalPositions();
	}

}

void Anchor::addTranslationParameters(RefinementStrategyPtr strategy,
                                      double mult)
{
	if (_motions.size() == 0)
	{
		return;
	}

	_motions[0]->addTranslationParameters(strategy);
}

void Anchor::addLibrationParameters(RefinementStrategyPtr strategy,
                                      int num)
{
	if (num < 0)
	{
		for (int i = 0; i < _quats.size(); i++)
		{
			_quats[i]->addVecToStrategy(strategy, 0.01, 0.0001, "rot");
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

		_quats[num]->addVecToStrategy(strategy, 0.01, 0.0001, "rot");
	}
}

void Anchor::addScrewParameters(RefinementStrategyPtr strategy,
                                int num)
{
	if (num < 0)
	{
		for (int i = 0; i < _quats.size(); i++)
		{
			_screws[i]->addVec2ToStrategy(strategy, 3.0, 0.01, "screw");
		}
	}
	else
	{
		if (num >= _quats.size())
		{
			return;
		}

		_screws[num]->addVec2ToStrategy(strategy, 3.0, 0.01, "screw");
	}
}

void Anchor::deleteLastScrew()
{
	if (_screws.size() == 0)
	{
		return;
	}

	_screws.pop_back();
	_quats.pop_back();
}
