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
	
	double rwork, rfree, swork, sfree;
	ave->averageRs(&rwork, &rfree, &swork, &sfree);
	
	VagFFTPtr fft = ave->getAverageFFT();
	
	std::vector<double> uc = fft->getUnitCell();

	f << "INSERT INTO Clusters" << std::endl;
	f << "(analysis_time, folder_path, average_rwork, average_rfree, "\
	"average_a, average_b, average_c, "\
	"average_alpha, average_beta, average_gamma) " << std::endl;

	/* analysis time */
	f << "SELECT FROM_UNIXTIME(" << i_to_str(t) << "), ";
	/* folder path */
	f << "'" << path << "', ";
	/* average Rwork */
	f << f_to_str(rwork, 3) + ", ";
	/* average Rfree */
	f << f_to_str(rfree, 3) + ", ";
	/* average uc-a */
	f << f_to_str(uc[0], 3) + ", ";
	/* average uc-b */
	f << f_to_str(uc[1], 3) + ", ";
	/* average uc-c */
	f << f_to_str(uc[2], 3) + ", ";
	/* average uc-alpha */
	f << f_to_str(uc[3], 3) + ", ";
	/* average uc-beta */
	f << f_to_str(uc[4], 3) + ", ";
	/* average uc-gamma */
	f << f_to_str(uc[5], 3) << std::endl;

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
