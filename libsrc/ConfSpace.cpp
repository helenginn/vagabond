// vagabond
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

#include "ConfSpace.h"
#include "SVDBond.h"
#include "ConfAxis.h"
#include "Polymer.h"
#include <hcsrc/FileReader.h>

ConfSpace::ConfSpace(int n)
{
	for (size_t i = 0; i < n; i++)
	{
		ConfAxis *axis = new ConfAxis();
		_axes.push_back(axis);
	}
}

void ConfSpace::calculateFrom(PolymerPtr pol)
{
	AtomList cas = pol->findAtoms("CA");

	_svd = new SVDBond(cas);
	_svd->setDoTorsion(true);
	_svd->performSVD();
	
	for (size_t i = 0; i < axisCount() && i < 16; i++)
	{
		axis(i)->getAxis(_svd, i);
	}

	assembleContributions();
}

void ConfSpace::readFromFile(std::string filename)
{
	std::string inside = get_file_contents(filename);
	
	std::vector<std::string> lines = split(inside, '\n');
	std::cout << lines.size() << std::endl;
	
	for (size_t i = 0; i < lines.size(); i++)
	{
		std::vector<std::string> bits = split(lines[i], ',');
		if (bits.size() == 0)
		{
			continue;
		}

		std::string atom = bits[0];
		std::vector<std::string> labels = split(atom, '_');
		
		if (labels.size() < 3)
		{
			continue;
		}
		
		int resi = atoi(labels[1].c_str());
		bool psi = (labels[2] == "psi");

		for (size_t j = 0; j < axisCount() && j < bits.size() - 1; j++)
		{
			std::string str = bits[j + 1];
			double val = atof(str.c_str());
			
			if (!psi)
			{
				axis(j)->setPhiDeviation(resi, val);
			}
			else
			{
				axis(j)->setPsiDeviation(resi, val);
			}
		}
	}

	assembleContributions();
}

void ConfSpace::assembleContributions()
{
	if (axisCount() == 0)
	{
		return;
	}

	if (axis(0)->count() == 0)
	{
		std::cout << "no residues" << std::endl;
		return;
	}
	
	int start = axis(0)->residueBegin();
	int end = axis(0)->residueEnd();

	for (size_t i = start; i <= end; i++)
	{
		std::vector<double> phis, psis;
		for (size_t j = 0; j < axisCount(); j++)
		{
			phis.push_back(axis(j)->getPhiDeviationForResidue(i));
			psis.push_back(axis(j)->getPsiDeviationForResidue(i));
		}
		
		_phis[i] = phis;
		_psis[i] = psis;
	}

}

ConfSpace::~ConfSpace()
{
	for (size_t i = 0; i < axisCount(); i++)
	{
		delete _axes[i];
	}
}

void ConfSpace::addMolecule(MoleculePtr mol)
{
	_molecules.push_back(mol);
}

double ConfSpace::totalMotionForResidue(int resi)
{
	double sum = 0;
	for (size_t i = 0; i < axisCount(); i++)
	{
		double a = axis(i)->getPsiDeviationForResidue(resi);
		a *= a;
		double b = axis(i)->getPhiDeviationForResidue(resi);
		b *= b;
		
		sum += a + b;
	}

	return sum;
}

double ConfSpace::weightedTorsion(int resi, bool phi, SpacePoint &sp)
{
	double sum = 0;
	std::vector<double> &torsions = (phi ? _phis[resi] : _psis[resi]);

	for (size_t i = 0; i < sp.size() && i < torsions.size(); i++)
	{
		sum += sp[i] * torsions[i];
	}

	return sum;
}
