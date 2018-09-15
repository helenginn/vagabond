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

Anchor::Anchor()
{
	_bFactor = 0;
	_absolute = empty_vec3();
}

std::vector<BondSample> *Anchor::getManyPositions()
{
	std::vector<BondSample> *bondSamples = &_storedSamples;
	bondSamples->clear();

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

	for (size_t i = 0; i < points.size(); i++)
	{
		vec3 full = vec3_add_vec3(points[i], _absolute);
		double occ = 1;
		occTotal += occ;

		BondSample sample;
		sample.basis = make_mat3x3();
		sample.occupancy = occ;
		sample.torsion = 0;
		sample.old_start = make_vec3(0, 0, 0);
		sample.start = full;

		bondSamples->push_back(sample);
	}

	for (size_t i = 0; i < bondSamples->size(); i++)
	{
		if (_occupancies.size() == bondSamples->size())
		{
			(*bondSamples)[i].occupancy = _occupancies[i];
		}
		else
		{
			(*bondSamples)[i].occupancy /= occTotal;
		}
	}

	return bondSamples;
}

void Anchor::addProperties()
{
	addDoubleProperty("bfactor", &_bFactor);
	addVec3Property("position", &_absolute);
	addReference("atom", _atom.lock());
	Model::addProperties();
}

std::string Anchor::getParserIdentifier()
{
	return "anchor_" + i_to_str(getAtom()->getAtomNum());
}

void Anchor::linkReference(ParserPtr object, std::string category)
{
	if (category == "atom")
	{
		AtomPtr atom = ToAtomPtr(object);
		_atom = atom;
	}
}
