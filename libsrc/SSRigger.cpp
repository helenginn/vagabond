//
//  SSRigger.cpp
//  vagabond
//
//  Created by Helen Ginn on 18/07/2017.
//  Copyright (c) 2018 Strubi. All rights reserved.
//

#include "SSRigger.h"
#include "Shouter.h"
#include "Crystal.h"
#include "Molecule.h"
#include "Monomer.h"
#include "Atom.h"

SSRigger::SSRigger()
{

}

void SSRigger::findCysteineSulphurs()
{
	for (size_t i = 0; i < _crystal->moleculeCount(); i++)
	{
		MoleculePtr molecule = _crystal->molecule(i);
		AtomList sulphurs = molecule->findAtoms("SG");
		
		/* Append to growing list */
		_cysSGs.reserve(_cysSGs.size() + sulphurs.size());

		for (size_t j = 0; j < sulphurs.size(); j++)
		{
			_cysSGs.push_back(sulphurs[j]);		
		}
	}
}

void SSRigger::convertCysteine(AtomPtr oneAtom)
{
	MonomerPtr monomer = oneAtom->getMonomer();

	for (size_t i = 0; i < monomer->atomCount(); i++)
	{
		AtomPtr atom = monomer->atom(i);
		atom->convertToDisulphide();	
	}
	
	monomer->propagateChange();
}

void SSRigger::findCloseCysteines()
{
	for (size_t i = 0; i < _cysSGs.size() - 1; i++)
	{
		for (size_t j = i + 1; j < _cysSGs.size(); j++)
		{
			AtomPtr aSulphur = _cysSGs[i];
			AtomPtr bSulphur = _cysSGs[j];	
			
			/* Atoms may be of different conformers of the same residue */
			if (aSulphur->getResidueNum() == bSulphur->getResidueNum())
			{
				std::cout << "\t" << aSulphur->shortDesc() << " and ";
				std::cout << bSulphur->shortDesc(); 
				std::cout << " are from the same residue." << std::endl;
				continue;	
			}
			
			/* Otherwise, if 2.6 Ang or less, consider disulphide. */
			bool close = aSulphur->closeToAtom(bSulphur, 2.6);
			
			if (close)
			{
				std::cout << "\tDetected disulphide between ";
				std::cout << aSulphur->shortDesc() << " and ";
				std::cout << bSulphur->shortDesc() << "." << std::endl;
				
				convertCysteine(aSulphur);
				convertCysteine(bSulphur);
				
				/* we do not remove them from the array */
			}
			else
			{
				std::cout << "\tNo disulphide between ";
				std::cout << aSulphur->shortDesc() << " and ";
				std::cout << bSulphur->shortDesc() << "." << std::endl;
			}
		}	
	}	
}

void SSRigger::findDisulphides()
{
	if (!_crystal)
	{
		shout_at_helen("No crystal before finding disulphides");
	}

	findCysteineSulphurs();
	
	std::cout << "\tCysteine sulphurs: " << _cysSGs.size() << std::endl;

	findCloseCysteines();
}
