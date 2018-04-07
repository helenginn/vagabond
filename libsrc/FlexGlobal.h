//
//  FlexGlobal.hpp
//  vagabond
//
//  Created by Helen Ginn on 27/12/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#ifndef FlexGlobal_hpp
#define FlexGlobal_hpp

#include <stdio.h>
#include "shared_ptrs.h"

typedef enum
{
	FlexTargetMaximiseIsotropy,
	FlexTargetMatchOrigBFactor,
	FlexTargetMatchElectronDensity,
} FlexTarget;

/**
 * \class FlexGlobal
 * \brief Looks after scoring function for whole molecule movements.
 */

class FlexGlobal
{
public:
	FlexGlobal();

	/** Set atoms to consider in the target function. */
	void setAtomGroup(AtomGroupPtr group)
	{
		_atomGroup = group;
	}

	static double score(void *object);

	/** Try to remove all anisotropic components, also use
	* setTargetBFactor(double) to choose the target B factor. */
	void maximiseIsotropy();

	/** Match the electron density of the AtomGroup using correlation
	* 	as a target function. */
	void matchElectronDensity()
	{
		_targetType = FlexTargetMatchElectronDensity;
	}

	/** Tries to match the anisotropic tensor to that found in the PDB file.
	* [Not recommended]
	*/
	void matchOriginalBees()
	{
		_targetType = FlexTargetMatchOrigBFactor;
	}

	/** Sets the target B factor if also using maximiseIsotropy(). */
	void setTargetBFactor(double value)
	{
		_targetIsoB = value;
	}
	
	/** Crystal against which electron density should be matched. */
	void setCrystal(CrystalPtr crystal)
	{
		_crystal = crystal;	
	}
private:
	double _targetIsoB;

	double matchOriginalBeeScore();
	double maximiseIsotropyScore();
	double matchElectronDensityScore();
	AtomGroupPtr _atomGroup;
	CrystalPtr _crystal;

	FlexTarget _targetType;
};

#endif /* FlexGlobal_hpp */
