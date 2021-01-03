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

std::string identifier(AtomPtr earlier, AtomPtr major, AtomPtr minor)
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
	
	std::string m = " ";
	if (minor)
	{
		m = minor->getAtomName();
	}

	return e + M + m;
}

AngleType Angler::mainAngleType(AtomPtr earlier, AtomPtr major, AtomPtr minor)
{
	AngleType type = AngleNone;
	std::string id = identifier(earlier, major, minor);

	if (false && _phiBond->getMajor()->getResidueNum() == 468)
	{
		std::cout << id << std::endl;
	}

	if (id == "CCAN" || id == "NCAC")
	{
		type = AngleNAC; 
	}
	else if (id == "CNCA" || id == "CANC")
	{
		type = AngleCNA;
	}
	else if (id == "CACN" || id == "NCCA")
	{
		type = AngleACN;
	}
	else if (id == "NCACB" || id == "CBCAN")
	{
		type = AngleNAB;
	}
	else if (id == "CCACB" || id == "CBCAC")
	{
		type = AngleBAC;
	}
	else if (id == "NCO" || id == "OCN")
	{
		type = AngleOCN;
	}
	else if (id == "CACO" || id == "OCCA")
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
		if (res == ResiNonPGIV)
		{
			if (ang == AngleNAC)
			{
				*where = GeomVariant::pointerToNonPGIVXProNAC();
			}
			if (ang == AngleCNA)
			{
				*where = GeomVariant::pointerToNonPGIVXProCNA();
			}
			if (ang == AngleNAB)
			{
				*where = GeomVariant::pointerToNonPGIVXProNAB();
			}
			if (ang == AngleBAC)
			{
				*where = GeomVariant::pointerToNonPGIVXProBAC();
			}
			if (ang == AngleOCN)
			{
				*where = GeomVariant::pointerToNonPGIVXProOCN();
			}
			if (ang == AngleACO)
			{
				*where = GeomVariant::pointerToNonPGIVXProACO();
			}
		}
		else if (res == ResiIleVal)
		{
			if (ang == AngleNAC)
			{
				*where = GeomVariant::pointerToIleValXProNAC();
			}
			if (ang == AngleCNA)
			{
				*where = GeomVariant::pointerToIleValXProCNA();
			}
			if (ang == AngleNAB)
			{
				*where = GeomVariant::pointerToIleValXProNAB();
			}
			if (ang == AngleBAC)
			{
				*where = GeomVariant::pointerToIleValXProBAC();
			}
			if (ang == AngleOCN)
			{
				*where = GeomVariant::pointerToIleValXProOCN();
			}
			if (ang == AngleACO)
			{
				*where = GeomVariant::pointerToIleValXProACO();
			}
		}
		else if (res == ResiGlycine)
		{
			if (ang == AngleNAC)
			{
				*where = GeomVariant::pointerToGlyXProNAC();
			}
			if (ang == AngleCNA)
			{
				*where = GeomVariant::pointerToGlyXProCNA();
			}
			if (ang == AngleNAB)
			{
				*where = GeomVariant::pointerToGlyXProNAB();
			}
			if (ang == AngleBAC)
			{
				*where = GeomVariant::pointerToGlyXProBAC();
			}
			if (ang == AngleOCN)
			{
				*where = GeomVariant::pointerToGlyXProOCN();
			}
			if (ang == AngleACO)
			{
				*where = GeomVariant::pointerToGlyXProACO();
			}
		}
		else if (res == ResiProline)
		{
			if (ang == AngleNAC)
			{
				*where = GeomVariant::pointerToProXProNAC();
			}
			if (ang == AngleCNA)
			{
				*where = GeomVariant::pointerToProXProCNA();
			}
			if (ang == AngleNAB)
			{
				*where = GeomVariant::pointerToProXProNAB();
			}
			if (ang == AngleBAC)
			{
				*where = GeomVariant::pointerToProXProBAC();
			}
			if (ang == AngleOCN)
			{
				*where = GeomVariant::pointerToProXProOCN();
			}
			if (ang == AngleACO)
			{
				*where = GeomVariant::pointerToProXProACO();
			}
		}
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
		else if (res == ResiGlycine)
		{
			if (ang == AngleNAC)
			{
				*where = GeomVariant::pointerToGlyNonXProNAC();
			}
			if (ang == AngleCNA)
			{
				*where = GeomVariant::pointerToGlyNonXProCNA();
			}
			if (ang == AngleNAB)
			{
				*where = GeomVariant::pointerToGlyNonXProNAB();
			}
			if (ang == AngleBAC)
			{
				*where = GeomVariant::pointerToGlyNonXProBAC();
			}
			if (ang == AngleOCN)
			{
				*where = GeomVariant::pointerToGlyNonXProOCN();
			}
			if (ang == AngleACO)
			{
				*where = GeomVariant::pointerToGlyNonXProACO();
			}
		}
		else if (res == ResiProline)
		{
			if (ang == AngleNAC)
			{
				*where = GeomVariant::pointerToProNonXProNAC();
			}
			if (ang == AngleCNA)
			{
				*where = GeomVariant::pointerToProNonXProCNA();
			}
			if (ang == AngleNAB)
			{
				*where = GeomVariant::pointerToProNonXProNAB();
			}
			if (ang == AngleBAC)
			{
				*where = GeomVariant::pointerToProNonXProBAC();
			}
			if (ang == AngleOCN)
			{
				*where = GeomVariant::pointerToProNonXProOCN();
			}
			if (ang == AngleACO)
			{
				*where = GeomVariant::pointerToProNonXProACO();
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

	AtomPtr central = _phiBond->getMajor();
	std::string res = central->getMonomer()->getIdentifier();
	
	_resi = ResiNonPGIV;

	if (res == "gly")
	{
		_resi = ResiGlycine;
	}
	
	if (res == "pro")
	{
		_resi = ResiProline;
	}
	
	if (res == "ile" || res == "val")
	{
		_resi = ResiIleVal;
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
	
	AngleType main = mainAngleType(earlier, major, minor);
	assignTable(_resi, main, &_table);
	AngleType sis = mainAngleType(earlier, major, sister);
	assignTable(_resi, sis, &_sisterTable);
	AngleType other = mainAngleType(minor, major, sister);
	assignTable(_resi, other, &_otherTable);

	if (false && _phiBond->getMajor()->getResidueNum() == 468)
	{
		std::cout << "Sister: " << sister << std::endl;
		std::cout << "main is " << main << std::endl;
		std::cout << "sis is " << sis << std::endl;
		std::cout << "other is " << other << std::endl;
	}
	
	if (_table != NULL)
	{
		_angledBond->setAngler(shared_from_this());
		_phiBond->addIndirectAngler(shared_from_this());
		_psiBond->addIndirectAngler(shared_from_this());
	}
	
	return (_table != NULL);
}

std::string Angler::getParserIdentifier()
{
	std::string result = getClassName() + "_" + _angledBond->shortDesc();
	result += "_";
	
	switch (_resi)
	{
		case ResiNonPGIV:
		result += "npgiv";
		break;
		case ResiIleVal:
		result += "npgiv";
		break;
		case ResiProline:
		result += "npgiv";
		break;
		case ResiGlycine:
		result += "npgiv";
		break;
		default:
		break;
	}

	result += "_";
	result += (_nextIsPro ? "pro" : "xpro");
	return result;
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
		if (_phiBond->getMajor()->getResidueNum() == 468)
		{
//			if (_angledBond->getMinor()->getAtomName() == "CB")
			{
				std::cout << _angledBond->shortDesc() << " ";
				std::cout << tPhi << " " << tPsi << " " << 
				" chosen " << angle << std::endl;
				std::cout << getParserIdentifier() << std::endl;
			}
		}
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

void Angler::forcePropagation()
{
	_angledBond->propagateChange(2);
	_phiBond->propagateChange(2);
	_psiBond->propagateChange(2);
}
