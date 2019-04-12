//
//  Anchor.cpp
//  vagabond
//
//  Created by Helen Ginn 2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "RefineMat3x3.h"
#include "Anchor.h"
#include "Absolute.h"
#include "Anisotropicator.h"
#include "Options.h"
#include "Crystal.h"
#include "Whack.h"
#include "Twist.h"
#include "mat4x4.h"
#include <sstream>

Anchor::Anchor(AbsolutePtr absolute) : ExplicitModel()
{
	_bFactor = absolute->getBFactor();
	_position = absolute->getAbsolutePosition();
	_molecule = absolute->getMolecule();
	_atom = absolute->getAtom();
	_trans = RefineMat3x3Ptr(new RefineMat3x3(this, cleanup));
	_libration = RefineMat3x3Ptr(new RefineMat3x3(this, cleanup));
	_screw = RefineMat3x3Ptr(new RefineMat3x3(this, cleanup));
	_libration->setZero();
	_screw->setZero();
	_disableWhacks = false;
	_alpha = 0;
	_beta = 0;
	_gamma = 0;
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
	_bFactor = 0;
	_alpha = 0;
	_beta = 0;
	_gamma = 0;
	_position = empty_vec3();
	_nDir = empty_vec3();
	_cDir = empty_vec3();
	_trans = RefineMat3x3Ptr(new RefineMat3x3(this, cleanup));
	_libration = RefineMat3x3Ptr(new RefineMat3x3(this, cleanup));
	_libration->setZero();
	_screw = RefineMat3x3Ptr(new RefineMat3x3(this, cleanup));
	_screw->setZero();
	_disableWhacks = false;
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

void Anchor::createStartPositions(Atom *callAtom)
{
	_storedSamples.clear();
	
	/* B factor isotropic only atm, get mean square displacement in
	 * each dimension. */
	double meanSqDisp = getBFactor() / (8 * M_PI * M_PI);
	meanSqDisp = sqrt(meanSqDisp);

	double occTotal = 0;

	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	int totalPoints = crystal->getSampleNum();
	totalPoints -= 1;
	
	if (totalPoints < 0)
	{
		totalPoints = 0;
	}

	double totalSurfaces = 0;
	int layers = 10;
	
	if (totalPoints < 20)
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

	int rnd = 1;
	std::vector<vec3> points;
	double increment = M_PI * (3.0 - sqrt(5));

	_sphereAngles.clear();
	vec3 sum = empty_vec3();
	double count = 0;

	for (int j = 0; j < layers; j++)
	{
		double m = meanSqDisp * (double)(j + 1) / (double)layers;

		int samples = layerSurfaces[j] * scale + 1;
		double offset = 2. / (double)samples;

		for (int i = 0; i < samples; i++)
		{
			double y = (((double)i * offset) - 1) + (offset / 2);
			double r = sqrt(1 - y * y);

			double phi = (double)((i + rnd) % samples) * increment;

			double x = cos(phi) * r;
			double z = sin(phi) * r;

			vec3 point = make_vec3(x * m, y * m, z * m);
			vec3_add_to_vec3(&sum, point);

			points.push_back(point);
			_sphereAngles.push_back(point);
			count++;
		}
	}
	
	vec3_mult(&sum, 1 / count);
	
	for (int i = 0; i < _sphereAngles.size(); i++)
	{
		vec3_subtract_from_vec3(&_sphereAngles[i], sum);
	}
	
	bool isN = (callAtom == &*(_nAtom.lock()));
	
	/* Get the rotation matrix for alpha, beta, gamma modifications */
	mat3x3 rot = getAnchorRotation();

	/* Want the direction to be the opposite of the calling bond */
	vec3 direction = isN ? _cDir : _nDir;
	vec3 other = isN ? _cDir2 : _nDir2;
	
	mat3x3_mult_vec(rot, &direction);
	mat3x3_mult_vec(rot, &other);
	
	vec3 empty = empty_vec3();

	for (size_t i = 0; i < points.size(); i++)
	{
		vec3 full = vec3_add_vec3(points[i], _position);
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

void Anchor::rotateBases()
{
	vec3 empty = empty_vec3();
	mat3x3 libration = _libration->getMat3x3();
	Anisotropicator tropicator;
	tropicator.setTensor(libration);
	/* Get the orthogonal basis vectors of this tensor */
	mat3x3 lib = tropicator.basis();

	/* Screw in units: Angstroms times degree */
	mat3x3 screw = _screw->getMat3x3();
	
	for (int i = 0; i < _storedSamples.size(); i++)
	{
		vec3 start = _storedSamples[i].start;
		vec3 neg_start = vec3_mult(start, -1);
		vec3 diff = vec3_subtract_vec3(start, _position);
		mat3x3 rot_only = make_mat3x3();

		for (int j = 0; j < 3; j++)
		{
			vec3 rot_vec = mat3x3_axis(lib, j);
			double dot = vec3_dot_vec3(diff, rot_vec);
			mat3x3 rot_mat = mat3x3_unit_vec_rotation(rot_vec, dot);
			mat3x3 basis = mat3x3_mult_mat3x3(rot_mat, 
			                                  _storedSamples[i].basis); 
			
			rot_only = mat3x3_mult_mat3x3(rot_mat, rot_only);
			_storedSamples[i].basis = basis;
		}
		
		vec3 diff_to_old = vec3_subtract_vec3(_storedSamples[i].old_start,
		                                      start);
		mat3x3_mult_vec(rot_only, &diff_to_old);
		vec3 new_old = vec3_add_vec3(start, diff_to_old);
		_storedSamples[i].old_start = new_old;
	}
}

mat3x3 Anchor::getAnchorRotation()
{
	mat3x3 rot = mat3x3_rotate(_alpha, _beta, _gamma);
	
	return rot;
}

void Anchor::translateStartPositions()
{
	mat3x3 translation = _trans->getMat3x3();
	Anisotropicator tropicator;
	tropicator.setTensor(translation);
	mat3x3 trans = tropicator.basis();

	vec3 sum = empty_vec3();

	for (int i = 0; i < _storedSamples.size(); i++)
	{
		vec3 start = _storedSamples[i].start;
		vec3 diff = vec3_subtract_vec3(start, _position);
		vec3 moved_diff = mat3x3_mult_vec(trans, diff);
		vec3 diffdiff = vec3_subtract_vec3(moved_diff, diff);
		
		vec3_add_to_vec3(&sum, diffdiff);

		vec3_add_to_vec3(&_storedSamples[i].start, diffdiff);
		vec3_add_to_vec3(&_storedSamples[i].old_start, diffdiff);
	}
	
	vec3_mult(&sum, -1 / (double)_storedSamples.size());

	for (int i = 0; i < _storedSamples.size(); i++)
	{
		vec3_add_to_vec3(&_storedSamples[i].start, sum);
		vec3_add_to_vec3(&_storedSamples[i].old_start, sum);
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
		_whacks[i]->getBond()->calculateMagicMat();
	}

	propagateChange(-1, true);

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
	
	createStartPositions(callAtom);

	/* Check if number of samples has changed for any reason - if so,
	 * initiate re-caching of initial atom positions for each Whack. */

	for (int i = 0; i < _whacks.size() && !_disableWhacks; i++)
	{
		bool refresh = _whacks[i]->needsRefresh(_storedSamples);	

		if (refresh)
		{
			recalculateWhacks();
			break;
		}
	}
	
	rotateBases();
	translateStartPositions();

	fixCentroid();

	
	/* Apply whacks as normal, if we are not re-caching Whacks. */
	for (int i = 0; i < _whacks.size() && !_disableWhacks; i++)
	{
		WhackPtr whack = _whacks[i];
		whack->applyToAnchorSamples(_storedSamples);
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
	addMat3x3Property("translation", _trans->getMat3x3Ptr());
	addMat3x3Property("libration", _libration->getMat3x3Ptr());
	addMat3x3Property("screw", _screw->getMat3x3Ptr());
	addReference("c_atom", _cAtom.lock());
	
	for (int i = 0; i < whackCount(); i++)
	{
		addChild("whack", _whacks[i]);
	}
	
	Model::addProperties();
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

void Anchor::linkReference(ParserPtr object, std::string category)
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

}

void Anchor::setPolymerBasis(mat3x3 basis)
{
	_libMotion = mat3x3_inverse(basis);

	double sum = 0;
	for (int i = 0; i < 9; i++)
	{
		sum += _libMotion.vals[i];
	}
	
	sum /= 9;
	mat3x3_mult_scalar(&_libMotion, 0.10 / sum);
	
	for (int i = 0; i < 9; i++)
	{
		_libMotion.vals[i] = std::min(0.5, _libMotion.vals[i]);
	}
	
}

void Anchor::addTranslationParameters(RefinementStrategyPtr strategy,
                                      double mult)
{
	_trans->addTensorToStrategy(strategy, 0.5 * mult, 0.001, "tr");
}

void Anchor::addLibrationParameters(RefinementStrategyPtr strategy,
                                      double mult)
{
	_libration->addTensorToStrategy(strategy, 0.2 * mult, 0.005, "li");
}

void Anchor::addScrewParameters(RefinementStrategyPtr strategy,
                                double mult)
{
	_screw->addMatrixToStrategy(strategy, 0.2 * mult, 0.005, "sc");
}

