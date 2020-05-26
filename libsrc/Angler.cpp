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

#include "Angler.h"
#include "Monomer.h"
#include "Atom.h"
#include "Bond.h"
#include "GeomVariant.h"
#include "Options.h"

Angler::Angler()
{
	_nextIsPro = false;
	_table = NULL;
	_sisterTable = NULL;
	_otherTable = NULL;
}

std::string identifier(AtomPtr earlier, AtomPtr major,
                       AtomPtr sister, AtomPtr minor)
{
	std::string e = " ";
	if (earlier)
	{
		e = earlier->getAtomName();
	}

	std::string M = " ";
	if (major)
	{
		M = major->getAtomName();
	}

	std::string s = " ";
	if (sister)
	{
		s = sister->getAtomName();
	}
	
	std::string m = " ";
	if (minor)
	{
		m = minor->getAtomName();
	}

	return e + M + m + s;
}

AngleType Angler::mainAngleType(AtomPtr earlier, AtomPtr major,
                                AtomPtr sister, AtomPtr minor)
{
	AngleType type = AngleNone;
	std::string id = identifier(earlier, major, AtomPtr(), minor);

	if (id == "CCAN " || id == "NCAC ")
	{
		type = AngleNAC; 
	}
	else if (id == "CNCA " || id == "CANC ")
	{
		type = AngleCNA;
	}
	else if (id == "CACN " || id == "NCCA ")
	{
		type = AngleACN;
	}
	else if (id == "NCACB " || id == "CBCAN ")
	{
		type = AngleNAB;
	}
	else if (id == "CCACB " || id == "CBCAC ")
	{
		type = AngleBAC;
	}
	else if (id == "NCO " || id == "OCN ")
	{
		type = AngleOCN;
	}
	else if (id == "CACO " || id == "OCCA ")
	{
		type = AngleACO;
	}

	return type;
}

void Angler::assignTable(ResiType res, AngleType ang, double **where)
{
	if (ang == AngleNone)
	{
		return;
	}

	if (_nextIsPro)
	{

	}
	else
	{
		if (res == ResiNonPGIV)
		{
			if (ang == AngleNAC)
			{
				*where = GeomVariant::pointerToNonPGIVNonXProNAC();
			}
			if (ang == AngleCNA)
			{
				*where = GeomVariant::pointerToNonPGIVNonXProCNA();
			}
			if (ang == AngleNAB)
			{
				*where = GeomVariant::pointerToNonPGIVNonXProNAB();
			}
			if (ang == AngleBAC)
			{
				*where = GeomVariant::pointerToNonPGIVNonXProBAC();
			}
			if (ang == AngleOCN)
			{
				*where = GeomVariant::pointerToNonPGIVNonXProOCN();
			}
			if (ang == AngleACO)
			{
				*where = GeomVariant::pointerToNonPGIVNonXProACO();
			}
		}
		else if (res == ResiIleVal)
		{
			if (ang == AngleNAC)
			{
				*where = GeomVariant::pointerToIleValNonXProNAC();
			}
			if (ang == AngleCNA)
			{
				*where = GeomVariant::pointerToIleValNonXProCNA();
			}
			if (ang == AngleNAB)
			{
				*where = GeomVariant::pointerToIleValNonXProNAB();
			}
			if (ang == AngleBAC)
			{
				*where = GeomVariant::pointerToIleValNonXProBAC();
			}
			if (ang == AngleOCN)
			{
				*where = GeomVariant::pointerToIleValNonXProOCN();
			}
			if (ang == AngleACO)
			{
				*where = GeomVariant::pointerToIleValNonXProACO();
			}
		}
	}
}

bool Angler::setupTable()
{
	if (!_angledBond || !_phiBond || !_psiBond)
	{
		return false;
	}

	if (!Options::isAngling())
	{
		return false;
	}

	AtomPtr major = _angledBond->getMajor();
	AtomPtr minor = _angledBond->getMinor();
	
	if (!major || !minor || !minor->getMonomer())
	{
		return false;
	}

	std::string res = minor->getMonomer()->getIdentifier();
	
	ResiType resi = ResiNonPGIV;

	if (res == "gly" || res == "pro")
	{
		return false;
	}
	
	if (res == "ile" || res == "val")
	{
		resi = ResiIleVal;
	}

	ModelPtr upstream = _angledBond->getParentModel();
	
	if (!upstream->isBond())
	{
		return false;
	}

	AtomPtr earlier = ToBondPtr(upstream)->getMajor();
	BondPtr sisterBond = _angledBond->getClosestSister();
	AtomPtr sister;
	
	if (sisterBond)
	{
		sister = sisterBond->getMinor();
	}
	
	AngleType main = mainAngleType(earlier, major, sister, minor);
	assignTable(resi, main, &_table);
	AngleType sis = mainAngleType(minor, major, sister, minor);
	assignTable(resi, sis, &_sisterTable);
	AngleType other = mainAngleType(earlier, major, minor, sister);
	assignTable(resi, other, &_otherTable);
	
	if (_table != NULL)
	{
		_angledBond->setAngler(shared_from_this());
	}
	
	return (_table != NULL);
}

std::string Angler::getParserIdentifier()
{
	return getClassName() + "_" + _angledBond->shortDesc() + 
	"_" + (_nextIsPro ? "_pro" : "_xpro");
}

void Angler::addProperties()
{

}

void Angler::addObject(ParserPtr object, std::string category)
{

}

void Angler::postParseTidy()
{

}

double Angler::getAngle(bool report, double *which)
{
	if (which == NULL)
	{
		return -1;
	}

	double tPhi = rad2deg(Bond::getTorsion(&*_phiBond));
	double tPsi = rad2deg(Bond::getTorsion(&*_psiBond));
	
	while (tPhi >= 180)
	{
		tPhi -= 360;
	}

	while (tPhi < -180)
	{
		tPhi += 360;
	}
	
	while (tPsi >= 180)
	{
		tPsi -= 360;
	}

	while (tPsi < -180)
	{
		tPsi += 360;
	}
	
	double angle = GeomVariant::interpolateAngle(which, tPhi, tPsi);
	
	if (report)
	{
		std::cout << tPhi << " " << tPsi << " " << angle << std::endl;
	}

	return deg2rad(angle);
}

double Angler::getMainAngle()
{
	return getAngle(false, _table);
}

double Angler::getSisterAngle()
{
	return getAngle(false, _sisterTable);
}

double Angler::getOtherAngle()
{
	return getAngle(false, _otherTable);
}
