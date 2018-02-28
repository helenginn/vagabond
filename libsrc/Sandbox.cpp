//
//  Sandbox.cpp
//  vagabond
//
//  Created by Helen Ginn on 18/08/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Sandbox.h"
#include <math.h>
#include "shared_ptrs.h"
#include <map>
#include "CSV.h"

void outputTorsionStuff()
{
	CSVPtr csv = CSVPtr(new CSV(3, "r1", "r2", "val"));

	double t1 = deg2rad(0);
	double t2 = deg2rad(0);
	double r1, r2;

	double cost1 = cos(t1);
	double cost2 = cos(t2);
	double sint1 = sin(t1);
	double sint2 = sin(t2);

	for (r1 = 0; r1 < deg2rad(360); r1 += deg2rad(2))
	{
		for (r2 = 0; r2 < deg2rad(360); r2 += deg2rad(2))
		{
			double cosr1 = cos(r1);
			double cosr2 = cos(r2);
			double sinr1 = sin(r1);
			double sinr2 = sin(r2);

			double distSq = 0;
			distSq += pow(cost1 - cosr1, 2);
			distSq += pow(1 - cost1 - cost2 - (1 - cosr1 - cosr2), 2);
			distSq += pow(2 * sint1 + sint2 - (2 * sinr1 + sinr2), 2);

			csv->addEntry(3, rad2deg(r1), rad2deg(r2), sqrt(distSq));
		}
	}

	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = "test";
	plotMap["height"] = "800";
	plotMap["width"] = "800";
	plotMap["xHeader0"] = "r1";
	plotMap["yHeader0"] = "r2";
	plotMap["zHeader0"] = "val";

	plotMap["xTitle0"] = "r1 angle";
	plotMap["yTitle0"] = "r2 angle";
	plotMap["style0"] = "heatmap";
	plotMap["colour0"] = "blue";
	plotMap["stride"] = "180";

	csv->plotPNG(plotMap);
}
