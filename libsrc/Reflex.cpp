//
//  Reflex.cpp
//  vagabond
//
//  Created by Helen Ginn, 2018
//  Copyright © 2018 Helen Ginn. All rights reserved.
//

#include "Reflex.h"
#include "Shouter.h"
#include "Polymer.h"
#include "FileReader.h"
#include "Monomer.h"
#include "CSV.h"
#include "Backbone.h"
#include "Options.h"
#include "AtomGroup.h"
#include "Atom.h"
#include "RefinementNelderMead.h"

Reflex::Reflex()
{
	_pieceCount = 1;
	_pieceType = PieceTypeResidue;
	_bFactor = 0;
	_kAbs = 1;
	_intercept = 0;
	_gradient = 1;
	_aniso = RefineMat3x3Ptr(new RefineMat3x3(this, NULL));
	_aniso->setZero();
}

void Reflex::segmentPolymer()
{
	/* Clear atom target B factors */
	
	for (int i = 0; i < _polymer->atomCount(); i++)
	{
		AtomPtr atom = _polymer->atom(i);
		atom->clearTargetB();
	}
	
	/* Divide the polymer into sections */
	
	for (size_t i = _polymer->monomerBegin(); 
	     i <= _polymer->monomerEnd() - _pieceCount; i++)
	{
		AtomGroupPtr segment = AtomGroupPtr(new AtomGroup());

		for (int j = i; j < i + _pieceCount; j++)
		{
			MonomerPtr monomer = _polymer->getMonomer(j);
			
			if (!monomer)
			{
				continue;
			}
			
			segment->addAtomsFrom(monomer->getBackbone());
		}

		_segments.push_back(segment);
	}
}

double Reflex::bFactorScore(void *object)
{
	Reflex *me = static_cast<Reflex *>(object);
	double b = me->_bFactor;
	double abs = me->_kAbs;
	FFTPtr obs = me->_obs;
	FFTPtr calc = me->_calc;

	mat3x3 inv = obs->getBasis();	
	mat3x3_scale(&inv, obs->nx, obs->ny, obs->nz);
	mat3x3 basis = mat3x3_inverse(inv);
	DiffractionPtr data = Options::getRuntimeOptions()->getActiveData();
	double res = 0; // in Angstroms
	res = me->_workspace.crystal->getMaxResolution(data);
	res /= 2;

	std::vector<double> xs, ys;
	mat3x3 aniso = me->_aniso->getMat3x3();
	
	for (int k = 0; k < obs->nz; k++)
	{
		for (int j = 0; j < obs->ny; j++)
		{
			for (int i = 0; i < obs->nx; i++)
			{
				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(basis, &ijk);
				double length = vec3_length(ijk); // recip
				
				double d = 1 / length; // real
				if (d < res)
				{
					continue;
				}

				double four_d_sq = (4 * d * d);
				double isoBMod = exp(-b / four_d_sq);

				long ele = obs->element(i, j, k);

				float fCalc = sqrt(calc->getIntensity(ele)) * isoBMod * 
				me->_kAbs;
				float fObs = sqrt(obs->getIntensity(ele));

				xs.push_back(fCalc);
				ys.push_back(fObs);
			}
		}
	}
	
	double rFac = r_factor(xs, ys);

	return rFac;
}

void Reflex::prepareSegments(AtomGroupPtr segment, bool preprocess)
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	CSym::CCP4SPG *spg = CSym::ccp4spg_load_by_standard_num(1);
	
	_workspace.scoreType = ScoreTypeCorrel;
	_workspace.crystal = crystal;
	_workspace.selectAtoms = segment->getAtoms();
	_workspace.segment = FFTPtr();
	_workspace.ave = empty_vec3();
	_workspace.basis = make_mat3x3();
	_workspace.flag = MapScoreFlagNone;
	
	double score = AtomGroup::scoreWithMapGeneral(&_workspace);
	FFTPtr calcSegment = FFTPtr(new FFT(*_workspace.segment));

	calcSegment->shiftToCentre();
	
	_workspace.scoreType = ScoreTypeCopyToSmaller;
	AtomGroup::scoreWithMapGeneral(&_workspace);
	FFTPtr obsSegment = _workspace.segment;
	_calcSize = calcSegment->averageAll();

	if (preprocess)
	{
		obsSegment->fft(1);
		calcSegment->fft(1);
		
		for (int i = 0; i < calcSegment->nn; i++)
		{
			calcSegment->data[i][0] *= _gradient;
			calcSegment->data[i][1] *= _gradient;

			obsSegment->data[i][0] *= _gradient;
			obsSegment->data[i][1] *= _gradient;
		}
		
		_calcSize = obsSegment->averageAll();
	}
	
	_obs = obsSegment;
	_calc = calcSegment;
}

void Reflex::refineFlex()
{	
	_bFactor = 0;
	_kAbs = 1;
	_aniso->setZero();
	
	NelderMeadPtr nm;
	nm = NelderMeadPtr(new RefinementNelderMead());
	nm->setJobName("segment_" + i_to_str(_round));
	nm->setEvaluationFunction(bFactorScore, this);
	nm->setCycles(20);
	nm->setSilent();
	nm->addParameter(this, getBFactor, setBFactor, 5, 0.01, "b");
	nm->addParameter(this, getAbsoluteScale, setAbsoluteScale, 5, 0.01, "k");
	nm->refine();
	nm->clearParameters();
}

void Reflex::calculate()
{
	if (_pieceType != PieceTypeResidue)
	{
		shout_at_helen("Only support residue flex estimation");
	}
	
	segmentPolymer();
	
	CSVPtr csv = CSVPtr(new CSV(2, "seg", "calc_ave"));
	
	std::vector<double> refineBs;
	
	for (int i = 0; i < _segments.size(); i++)
	{
		_round = i + _pieceCount / 2;
		prepareSegments(_segments[i]);
		refineFlex();
		
		for (int j = 0; j < _segments[i]->atomCount(); j++)
		{
			_segments[i]->atom(j)->setTargetB(_bFactor);
		}
		
		refineBs.push_back(_bFactor);
		csv->addEntry(2, (double)_round, _bFactor);
	}
		
	csv->writeToFile(_polymer->getChainID() + "_bfactor_targets.csv");
}
