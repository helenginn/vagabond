// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#ifndef __vagabond__VBondReader__
#define __vagabond__VBondReader__

#include "shared_ptrs.h"

/**
 * \class VBondReader
 * \brief Takes a filename of .vbond format and generates the equivalent
 * Vagabond-style model corresponding to it.
 */

class VBondReader
{
public:
	VBondReader();

	void setFilename(std::string filename)
	{
		_filename = filename;
	}

	CrystalPtr getCrystal();
private:
	std::string _filename;
};


#endif
