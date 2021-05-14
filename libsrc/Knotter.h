//
//  Knotter.h
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Knotter__
#define __vagabond__Knotter__

#include <stdio.h>
#include "shared_ptrs.h"

/**
 * \class Knotter
 * \brief Links atoms together using known (hard-coded) topology of
 * protein primary sequence */

class Knotter
{
public:
	Knotter();

	void setSidechain(SidechainPtr sidechain)
	{
		_sidechain = sidechain;
	}

	void setBackbone(BackbonePtr backbone)
	{
		_backbone = backbone;
	}

	void tie();
	void tieTowardsCTerminus();
	void tieTowardsNTerminus();
	void betaAngler(bool onRight);
private:
	void makeAngler(BondPtr phi, BondPtr psi, MonomerPtr mon,
                         std::string atomName);

	AnglerPtr setupAngler(MonomerPtr mon, int resi);
	SidechainPtr _sidechain;
	BackbonePtr _backbone;

	BondPtr tieBetaCarbon(AtomPtr torsionAtom);
	
	int _bondAngles;

	void makeAlanine();
	void makeGlycine();
	void makeCysteine();
	void makeValine();
	void makeSerine();
	void makeLysine();
	void makeLeucine();
	void makeIsoleucine();
	void makeThreonine();
	void makeHistidine();
	void makeTyrosine();
	void makeAspartate();
	void makePhenylalanine();
	void makeTryptophan();
	void makeMethionine();
	void makeSelenoMet();
	void makeProline();
	void makeArginine();
	void makeAsparagine();
	void makeGlutamine();
	void makeGlutamate();
};


#endif /* defined(__vagabond__Knotter__) */
