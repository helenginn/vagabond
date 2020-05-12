// cluster4x
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

#include "Output.h"
#include <fstream>
#include <time.h>
#include <libsrc/FileReader.h>
#include "MtzFile.h"
#include "Group.h"

void Output::createSQL(Group *ave, std::string path)
{
	std::ofstream f;
	f.open(path + "/insert.sql");
	
	f << "-- SQL input for cluster" << std::endl;
	f << "-- All queries should be safe to execute more than once" << std::endl;
	f << std::endl;
	f << "-- set up initial Cluster entry if this does not exist" << std::endl;
	f << std::endl;
	
	time_t t;
	time(&t);

	f << "INSERT INTO Clusters" << std::endl;
	f << "(analysis_time, folder_path)" << std::endl;
	f << "SELECT FROM_UNIXTIME(" << i_to_str(t) << "), ";
	f << "'" << path << "'" << std::endl;
	f << "WHERE NOT EXISTS" << std::endl;
	f << "(SELECT 1 FROM Clusters WHERE folder_path = '" << path << "');";
	f << std::endl;
	f << std::endl;
	
	f << " -- linking refinement IDs to this cluster" << std::endl;
	f << std::endl;
	
	for (size_t i = 0; i < ave->mtzCount(); i++)
	{
		MtzFile *file = ave->getMtzFile(i);
		
		if (file->isDead())
		{
			continue;
		}
		
		std::string r_id = file->refinementID();

		f << "INSERT INTO Cluster_Members" << std::endl;
		f << "(cluster_id, refinement_id)" << std::endl;
		f << "SELECT cluster_id, " << r_id << " FROM Clusters WHERE ";
		f << "(folder_path = '" << path << "')" << std::endl;
		f << "AND NOT EXISTS" << std::endl;
		f << "(SELECT 1 FROM Clusters JOIN Cluster_Members ON" << std::endl;
		f << "Clusters.cluster_id = Cluster_Members.cluster_id" << std::endl;
		f << "WHERE folder_path = '" << path << "' AND refinement_id = ";
		f << r_id << ");" << std::endl;
		f << std::endl;
	}
	
	f << "-- done." << std::endl;

	f.close();

}
