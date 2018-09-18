//
//  Anchor.cpp
//  vagabond
//
//  Created by Helen Ginn 2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "Anchor.h"
#include "Absolute.h"
#include "Options.h"
#include "Crystal.h"

Anchor::Anchor(AbsolutePtr absolute)
{
	_bFactor = absolute->getBFactor();
	_absolute = absolute->getAbsolutePosition();
	_molecule = absolute->getMolecule();
	_atom = absolute->getAtom();
}

void Anchor::setNeighbouringAtoms(AtomPtr nPre, AtomPtr nAtom, 
                                  AtomPtr cAtom, AtomPtr cPost)
{
	_nAtom = nAtom;
	_cAtom = cAtom;
//	_nPre = nPre;
//	_cPost = cPost;
	
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
	_absolute = empty_vec3();
	_nDir = empty_vec3();
	_cDir = empty_vec3();
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
		vec3 full = vec3_add_vec3(points[i], _absolute);
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
	Atom *callAtom = static_cast<Atom *>(caller);
	createStartPositions(callAtom);
	return &_storedSamples;
}

void Anchor::addProperties()
{
	addDoubleProperty("bfactor", &_bFactor);
	addVec3Property("position", &_absolute);
	addReference("atom", _atom.lock());
	addReference("n_atom", _nAtom.lock());
	addVec3Property("n_dir", &_nDir);
	addVec3Property("c_dir", &_cDir);
	addVec3Property("pre_n", &_nDir2);
	addVec3Property("post_c", &_cDir2);
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
