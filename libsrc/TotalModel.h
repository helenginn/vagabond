// Vagabond
// Copyright (C) 2019 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#ifndef __Vagabond__TotalModel__
#define __Vagabond__TotalModel__

#include <stdio.h>
#include <vector>
#include "shared_ptrs.h"
#include <hcsrc/mat3x3.h>
#include <string>
#include <map>
#include <hcsrc/maths.h>
#include "AtomGroup.h"
#include "../libccp4/csymlib.h"
#include <iostream>

typedef std::map<std::string, MoleculePtr> MoleculeMap;

class ConfSpace;

class TotalModel : public AtomGroup
{
public:
	TotalModel();

	TotalModelPtr shared_from_this()
	{
		return ToTotalModelPtr(BaseParser::shared_from_this());
	}

	/**
	* How many molecules are included in a TotalModel.
	*/
	size_t moleculeCount()
	{
		return _molecules.size();
	}

	/**
	* 	get the stored molecule. Cannot guarantee order of molecules will be
	* 	the same as in the PDB file, but adequate for looping over molecules.
	*   \param i molecule number to obtain.
	* 	*/
	MoleculePtr molecule(long int i)
	{
		MoleculeMap::iterator it = _molecules.begin();
		std::advance(it, i);
		return it->second;
	}

	void addMolecule(MoleculePtr molecule);

	/**
	* 	get the stored molecule by chain name.
	*/
	MoleculePtr molecule(std::string chain)
	{
		if (_molecules.count(chain))
		{
			return _molecules[chain];
		}

		return MoleculePtr();
	}


	/** Calculates the anchor residue for each Polymer and assigns to each. */
	void setAnchors();
	
	size_t polymerCount();

	void refreshAnchors();
	void removeMolecule(MoleculePtr mol);
	void removeAtom(AtomPtr atom);

	void recalculateAtoms();
	
	void addMotion(MotionPtr mot, PolymerPtr pol = PolymerPtr());
	
	void setFilename(std::string file)
	{
		_filename = file;
	}
	
	std::string getFilename()
	{
		return _filename;
	}

	size_t motionCount()
	{
		return _motions.size();
	}

	void makeOverallMotion();
	MotionPtr getOverallMotion();
	void resetMotions();

	/** Obtain the current number of samples, i.e. number of conformers
	 * in the ensemble. */
	int getSampleNum();

	/** Total average B factor across model derived from Vagabond model */
	double averageBFactor();

	/** Prints scattering proportion in this crystal determined by bonds. */
	void tiedUpScattering();
	
	WaterNetworkPtr getWaterNetwork();
	void setupConformationalSpace();
protected:

	virtual void addProperties();
	virtual void postParseTidy();
	virtual void addObject(ParserPtr object, std::string category);

	MoleculeMap _molecules;
	std::vector<MotionPtr> _motions;
	std::string _filename;
	bool _tied;
	int _sampleNum;
	ConfSpace *_confSpace;
private:

};

#endif
