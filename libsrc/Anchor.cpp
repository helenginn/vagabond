//
//  Anchor.cpp
//  vagabond
//
//  Created by Helen Ginn 2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "Anchor.h"
#include "Absolute.h"
#include "Anisotropicator.h"
#include "Options.h"
#include "Crystal.h"
#include "mat4x4.h"
#include <sstream>

Anchor::Anchor(AbsolutePtr absolute)
{
	_bFactor = absolute->getBFactor();
	_position = absolute->getAbsolutePosition();
	_molecule = absolute->getMolecule();
	_atom = absolute->getAtom();
	_translation = make_mat3x3();
	_rotVec = make_vec3(0, 0, 0);
	_rotCentre = make_vec3(0, 0, 0);
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

Anchor::Anchor()
{
	_bFactor = 0;
	_position = empty_vec3();
	_nDir = empty_vec3();
	_cDir = empty_vec3();
	_rotVec = make_vec3(0, 0, 0);
	_rotCentre = make_vec3(0, 0, 0);
	_translation = make_mat3x3();
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

			points.push_back(point);
			_sphereAngles.push_back(point);
		}
	}
	
	bool isN = (callAtom == &*(_nAtom.lock()));

	/* Want the direction to be the opposite of the calling bond */
	vec3 *direction = isN ? &_cDir : &_nDir;
	vec3 *other = isN ? &_cDir2 : &_nDir2;
	vec3 empty = empty_vec3();

	for (size_t i = 0; i < points.size(); i++)
	{
		vec3 full = vec3_add_vec3(points[i], _position);
		vec3 next = vec3_add_vec3(full, *direction);
		vec3 prev = vec3_add_vec3(full, *other);
	
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
	if (vec3_length(_rotVec) < 1e-6)
	{
		return;
	}

	vec3 empty = empty_vec3();
	vec3 reverseCentre = _rotCentre;
	vec3_mult(&reverseCentre, -1);
	
	for (int i = 0; i < _storedSamples.size(); i++)
	{
		vec3 start = _storedSamples[i].start;
		vec3 diff = vec3_subtract_vec3(start, _position);
		double dot = vec3_dot_vec3(diff, _rotVec);
		mat3x3 rot_mat = mat3x3_unit_vec_rotation(_rotVec, dot);
		mat4x4 rot_mat4 = mat4x4_from_rot_trans(rot_mat, empty);

		mat4x4 change = make_mat4x4();
		mat4x4_translate(&change, _rotCentre);
		change = mat4x4_mult_mat4x4(rot_mat4, change);
		mat4x4_translate(&change, reverseCentre);

		mat4x4 basis = mat4x4_from_rot_trans(_storedSamples[i].basis, empty);
		mat4x4 transformed = mat4x4_mult_mat4x4(change, basis);

		mat3x3 tmp = mat4x4_get_rot(transformed);
		_storedSamples[i].basis = tmp;
	}
}

void Anchor::translateStartPositions()
{
	Anisotropicator tropicator;
	tropicator.setTensor(_translation);
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

std::vector<BondSample> *Anchor::getManyPositions(void *caller)
{
	std::vector<BondSample> copied = _storedSamples;
	Atom *callAtom = static_cast<Atom *>(caller);
	createStartPositions(callAtom);
	rotateBases();
	translateStartPositions();
	sanityCheck();
	
	return &_storedSamples;
}

void Anchor::addProperties()
{
	addDoubleProperty("bfactor", &_bFactor);
	addVec3Property("position", &_position);
	addReference("atom", _atom.lock());
	addReference("n_atom", _nAtom.lock());
	addVec3Property("n_dir", &_nDir);
	addVec3Property("c_dir", &_cDir);
	addVec3Property("pre_n", &_nDir2);
	addVec3Property("post_c", &_cDir2);
	addVec3Property("rot_vec", &_rotVec);
	addVec3Property("rot_centre", &_rotCentre);
	addMat3x3Property("translation", &_translation);
	addReference("c_atom", _cAtom.lock());
	Model::addProperties();
}

std::string Anchor::getParserIdentifier()
{
	return "anchor_" + i_to_str(getAtom()->getAtomNum());
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
	_cAtom.lock()->getModel()->propagateChange(depth, refresh);
	_nAtom.lock()->getModel()->propagateChange(depth, refresh);

	/* Will force recalculation of final positions */
	Model::propagateChange(depth, refresh);

}
