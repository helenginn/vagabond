//
//  DiffractionMTZ.h
//  vagabond
//
//  Created by Helen Ginn on 22/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__DiffractionMTZ__
#define __vagabond__DiffractionMTZ__

#include <stdio.h>
#include "Diffraction.h"
#include "Shouter.h"
#include "../libccp4/cmtzlib.h"

/**
 * \class DiffractionMtz
 * \brief Diffraction which obtains diffraction information from a MTZ file.
 */

class DiffractionMtz : public Diffraction
{
public:
	virtual void load();

	void syminfoCheck();
	
	void setDifference(std::string lp, std::string lm)
	{
		_plus = lp;
		_minus = lm;
	}
	
	void avoidPhases(bool t)
	{
		_avoidPhases = t;
	}
private:
	LabelChoice prepareChoice(CMtz::MTZ *mtz);

	bool _avoidPhases = false;
	std::string _plus;
	std::string _minus;
};

#endif /* defined(__vagabond__DiffractionMTZ__) */
