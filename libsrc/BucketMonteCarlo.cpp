//
//  BucketMonteCarlo.cpp
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "BucketMonteCarlo.h"
#include "BucketBulkSolvent.h"
#include "Crystal.h"
#include "CSV.h"
#include "mat3x3.h"
#include "Molecule.h"
#include "Plucker.h"
#include "Atom.h"
#include "Absolute.h"
#include "Options.h"

BucketMonteCarlo::BucketMonteCarlo()
{
	_loadedAnalysis = 0;
	_extra = AtomGroupPtr(new AtomGroup());
}

void pluckTwoAtoms(Plucker *pluck, Atom **one, Atom **two)
{
	*one = static_cast<Atom *>(pluck->pluck());
	*two = (*one)->pluckAnother();
}

void BucketMonteCarlo::loadAnalysis(std::string filename)
{
	if (_loadedAnalysis) return;
	std::cout << "Loading previous solvent analysis..." << std::endl;

	std::string contents = get_file_contents(filename);
	
	const int divide = 5;
	setupNodes(divide);

	std::vector<std::string> lines = split(contents, '\n');

	for (int i = 1; i < lines.size(); i++)
	{
		std::vector<std::string> components = split(lines[i], ',');
		
		if (components.size() < 4) continue;

		double cdoubles[4];
		for (int j = 0; j < components.size(); j++)
		{
			cdoubles[j] = atof(components[j].c_str());
		}
		
		_mdnode->addToNode(cdoubles[3], 3, cdoubles[0], cdoubles[1],
		                  cdoubles[2]);
	}

	makePluckers(6.);
	_loadedAnalysis = 1;
}

void BucketMonteCarlo::makePluckers(double distance)
{	
	const int split = 5;
	double step = (MAX_CHECK_DISTANCE - MIN_CHECK_DISTANCE) / pow(2, split);
	
	for (double v = MIN_CHECK_DISTANCE; v <= distance; v += step)
	{
		_mdnode->makePlucker(0, v, 0);
		Plucker *plucker = _mdnode->getPlucker(1);
		_pluckerMap[v] = plucker;
	}
}


void BucketMonteCarlo::addSolvent()
{
	/** Need the regular bulk solvent analysis to be done. */
	BucketPtr tmpBucket = BucketPtr(new BucketBulkSolvent());
	tmpBucket->setCrystal(getCrystal());
	tmpBucket->addSolvent();
	tmpBucket->processMaskedRegions();
	_maskedRegions = tmpBucket->getMaskedRegions();
	
	std::string file = Options::getSolventFile();
	loadAnalysis(file);

	CrystalPtr crystal = getCrystal();
	mat3x3 real2frac = crystal->getReal2Frac();

	_solvent = FFTPtr(new FFT(*crystal->getFFT()));
	_solvent->setAll(0);
	
	std::cout << "Making imaginary water networks." << std::endl;
	_baseGrp = AtomGroupPtr(new AtomGroup());
	
	for (int i = 0; i < crystal->moleculeCount(); i++)
	{
		if (crystal->molecule(i)->getClassName() != "Molecule")
		{
			continue;
		}
		
		_baseGrp->addAtomsFrom(crystal->molecule(i));
	}

	AtomGroupPtr all = allWaters();
	for (int i = 0; i < 5; i++)
	{
		Plucker *plucker = all->makePluckableWaters();
		makeNewWaters(plucker, 1000);
		all = allWaters();
		delete plucker;
	}

	vec3 offset = make_vec3(0, 0, 0);

	for (int i = 0; i < _extra->atomCount(); i++)
	{
		_extra->atom(i)->addToMap(_solvent, real2frac, offset);
	}
	
	FFTPtr crystFFT = crystal->getFFT();
	std::vector<double> obses, calcs;

	CSVPtr csv = CSVPtr(new CSV(2, "obs", "calc"));
	for (int i = 0; i < _solvent->nn; i++)
	{
		if (_maskedRegions->getMask(i) != 1)
		{
			continue;
		}

		double obs = crystFFT->data[i][0];
		double calc = _solvent->data[i][0];
		csv->addEntry(2, obs, calc);
		
		obses.push_back(obs);
		calcs.push_back(calc);
	}
	
	double correl = correlation(obses, calcs);
	std::cout << "Correlation: " << correl << std::endl;

	csv->writeToFile("cc_score.csv");

	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = "solv_real_space";
	plotMap["xHeader0"] = "obs";
	plotMap["yHeader0"] = "calc";
	plotMap["colour0"] = "black";

	plotMap["xTitle0"] = "Observed";
	plotMap["yTitle0"] = "Calculated";
	plotMap["style0"] = "scatter";
	csv->plotPNG(plotMap);
}

AtomGroupPtr BucketMonteCarlo::allWaters()
{
	AtomGroupPtr newGrp = AtomGroupPtr(new AtomGroup());
	newGrp->addAtomsFrom(_baseGrp);
	newGrp->addAtomsFrom(_extra);
	
	return newGrp;
}

AtomGroupPtr BucketMonteCarlo::makeNewWaters(Plucker *plucker, int total)
{
	CrystalPtr crystal = getCrystal();
	
	std::string path = FileReader::addOutputDirectory("waters.pdb");

	int count = 0;
	std::ofstream file;
	file.open(path);
	
	AtomGroupPtr additions = AtomGroupPtr(new AtomGroup());

	std::cout << "Adding waters" << std::flush;

	while (count < total)
	{
		Atom *one, *two;
		pluckTwoAtoms(plucker, &one, &two);
		
		double dist = one->getDistanceFrom(two);
		double right, angle;
		int success = getRandomValues(dist, &right, &angle);
//		success = bucket->getReallyRandomValues(dist, &right, &angle);
		
		if (!success) continue;
		
		AbsolutePtr absOne = ToAbsolutePtr(one->getModel());
		AbsolutePtr absTwo = ToAbsolutePtr(two->getModel());

		vec3 pos1 = absOne->getRandomPosition();
		vec3 pos2 = absTwo->getRandomPosition();
	
		vec3 diff = vec3_subtract_vec3(pos2, pos1);
		mat3x3 ortho = mat3x3_ortho_axes(diff);
	
		/* First model... random second axis angle */
		double twizzle = rand() / (double)RAND_MAX;
		mat3x3 rot = mat3x3_rotate(0, 0, deg2rad(twizzle));
		mat3x3 combo = mat3x3_mult_mat3x3(ortho, rot);
		
		double ratio = tan(deg2rad(angle) - M_PI / 2);
		vec3 protoPos = make_vec3(1, 0, ratio);
		vec3_set_length(&protoPos, right);
		mat3x3_mult_vec(combo, &protoPos);
		
		vec3 finalPos = vec3_add_vec3(protoPos, pos2);
		
		/*
		for (int i = 0; i < 3; i++)
		{
			double val = rand() / (double)RAND_MAX;
			val *= 100;
			*(&finalPos.x + i) = val;
		}
		*/
		
		if (finalPos.x != finalPos.x) continue;

		if (!isSolvent(finalPos))
		{
			continue;
		}
		
		double occ = one->getModel()->getEffectiveOccupancy() *
		             two->getModel()->getEffectiveOccupancy();
		
		occ *= 0.01;

		AbsolutePtr abs = AbsolutePtr(new Absolute(finalPos, 20, "O", occ));
		int num = additions->issueAtomNumber();
		abs->setIdentity(1, "Z", "HOH", "O", num);
		abs->setHeteroAtom(true);
		AtomPtr three = abs->makeAtom();
		additions->addAtom(three);
		
		if (!three) continue;
		
		count++;
		
		if (count % 10 == 0)
		{
			std::cout << "." << std::flush;
		}
	}
	
	_extra->addAtomsFrom(additions);
	
	file << _extra->getPDBContribution(PDBTypeAverage);

	std::cout << std::endl;
	
	file.close();
	
	/*
	std::cout << "Scoring against map..." << std::endl;
	
	double val = -_extra->scoreWithMap(ScoreTypeCorrel, crystal, true);
	std::cout << "Correlation with reference (extra only): " << val << std::endl;
	*/
	
	return _extra;

}
	
