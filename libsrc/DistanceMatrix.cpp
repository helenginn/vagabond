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

#include "DistanceMatrix.h"
#include "Atom.h"
#include "Polymer.h"
#include "Monomer.h"
#include "CSV.h"
#include "maths.h"

DistanceMatrix::DistanceMatrix(PolymerPtr pol)
{
	_pol = pol;

}

void DistanceMatrix::draw()
{
	CSVPtr csv = CSVPtr(new CSV(6, "resi", "resj", "dist", "x_dev",
	                            "y_dev", "z_dev"));
	
	for (int i = _pol->monomerBegin(); i < _pol->monomerEnd(); i++)
	{
		MonomerPtr moni = _pol->getMonomer(i);
		
		if (!moni)
		{
			continue;
		}

		AtomPtr cai = moni->findAtom("CA");
		
		if (!cai)
		{
			continue;
		}

		for (int j = _pol->monomerBegin(); j < _pol->monomerEnd(); j++)
		{
			MonomerPtr monj = _pol->getMonomer(j);

			if (!monj)
			{
				continue;
			}

			AtomPtr caj = monj->findAtom("CA");

			if (!caj)
			{
				continue;
			}

			vec3 diff;
			double dist = compareAtoms(cai, caj, &diff);
			csv->addEntry(5, (double)i, (double)j, dist,
			              diff.x, diff.y, diff.z);
		}
	}
	
	std::string filename = _pol->getChainID() + "_distance";

	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = filename;
	plotMap["height"] = "800";
	plotMap["width"] = "800";
	plotMap["xHeader0"] = "resi";
	plotMap["yHeader0"] = "resj";
	plotMap["zHeader0"] = "dist";

	plotMap["xTitle0"] = "resi";
	plotMap["yTitle0"] = "resj";
	plotMap["style0"] = "heatmap";
	plotMap["stride"] = i_to_str(_pol->monomerEnd() - _pol->monomerBegin());

	csv->writeToFile(filename + ".csv");
	csv->plotPNG(plotMap);
}

double DistanceMatrix::compareAtoms(AtomPtr a, AtomPtr b, vec3 *difference)
{
	std::vector<BondSample> aPos = a->getExplicitModel()->getFinalPositions();
	std::vector<BondSample> bPos = b->getExplicitModel()->getFinalPositions();
	
	vec3 aAve = a->getModel()->getAbsolutePosition();
	vec3 bAve = b->getModel()->getAbsolutePosition();
	
	CorrelData cd = empty_CD();
	CorrelData cd2 = empty_CD();
	
	for (int i = 0; i < aPos.size(); i++)
	{
		vec3 aDiff = vec3_subtract_vec3(aPos.at(i).start, aAve);
		vec3 bDiff = vec3_subtract_vec3(bPos.at(i).start, bAve);
		vec3 diff = vec3_subtract_vec3(aDiff, bDiff);
		double len = vec3_length(diff);
		add_to_CD(&cd, diff.x, diff.y);
		add_to_CD(&cd2, diff.z, len);
	}

	double xm, ym, zm, xs, ys, zs, ss, mm;
	means_stdevs_CD(cd, &xm, &ym, &xs, &ys);
	means_stdevs_CD(cd2, &zm, &mm, &zs, &ss);

	if (difference != NULL)
	{
		*difference = make_vec3(xs, ys, zs);
		std::cout << vec3_desc(*difference) << std::endl;
	}

	return ss;
}
