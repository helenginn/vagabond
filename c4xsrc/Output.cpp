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
#include "Group.h"
#include "AveDiffraction.h"
#include "MtzFile.h"
#include <libsrc/FileReader.h>
#include <unistd.h>
#include <iostream>
#include <fstream>

Output::Output()
{

}

bool Output::prepCluster(Group *ave)
{
	if (!ave->isMarked() || ave->isExported())
	{
		return false;
	}

	std::string folder = findNewFolder("cluster4x_");
	std::string cwd = getcwd(NULL, 0);
	std::string path = cwd + "/" + folder;

	DIR *dir = opendir(path.c_str());

	if (dir)
	{
		closedir(dir);
		std::cout << "Warning! " << folder << " already exists!" << std::endl;
	}
	else if (ENOENT == errno)
	{
		mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}

	for (size_t i = 0; i < ave->mtzCount(); i++)
	{
		MtzFile *file = ave->getMtzFile(i);
		
		if (file->isDead())
		{
			continue;
		}
		
		std::string metadata = file->metadata();

		std::string subpath = path + "/" + metadata;
		mkdir(subpath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		
		std::string pdbpath = subpath + "/final.pdb";
		std::string mtzpath = subpath + "/final.mtz";
		
		std::string mtzlink = file->getFilename();
		if (mtzlink[0] != '/')
		{
			mtzlink = cwd + "/" + file->getFilename();
		}
		
		std::string pdblink = file->getPdbPath();
		if (pdblink[0] != '/')
		{
			pdblink = cwd + "/" + file->getPdbPath();
		}
		
		symlink(pdblink.c_str(), pdbpath.c_str());
		symlink(mtzlink.c_str(), mtzpath.c_str());
		
		std::cout << "Created soft links for " << metadata << std::endl;
	}
	
	std::string pandda = path + "/pandda.sh";
	createPanDDAFile(pandda);
	
	createNotes(ave, path + "/notes.txt");
	
	ave->setExported(true);
	return true;
}

void Output::createNotes(Group *ave, std::string file)
{
	std::ofstream f;
	f.open(file);
	int dead = 0;
	int total = 0;
	CorrelData cd = empty_CD();

	for (size_t i = 0; i < ave->mtzCount(); i++)
	{
		MtzFile *file = ave->getMtzFile(i);
		
		if (file->isDead())
		{
			dead++;
			continue;
		}
		
		add_to_CD(&cd, file->getRWork(), file->getRFree());
		total++;
	}
	
	double mean_rwork, mean_rfree, stdev_rwork, stdev_rfree;

	means_stdevs_CD(cd, &mean_rwork, &mean_rfree, 
	                &stdev_rwork, &stdev_rfree);
	
	f << "Notes for cluster:" << std::endl;
	f << "\tNo. datasets: " << total << std::endl;
	f << "\tNo. discarded dead sets: " << dead << std::endl;
	f << "\tAverage Rwork: " << mean_rwork << " (s="
	<< stdev_rwork << ")" << std::endl;
	f << "\tAverage Rfree: " << mean_rfree << " (s="
	<< stdev_rfree << ")" << std::endl;
	
	std::string ucInfo = AveDiffraction::unitCellDesc(ave->getAverageFFT());

	f << "\tAverage unit cell: " << ucInfo << std::endl;

	f.close();
}

void Output::createPanDDAFile(std::string file)
{
	std::ofstream f;
	f.open(file);

	f << "#!/bin/bash" << std::endl;
	f << std::endl;

	f << "mycpus=`grep '^processor' /proc/cpuinfo | wc -l`" << std::endl;
	f << "cpus=`echo \"32 $mycpus\" | tr ' ' \"\n\" | sort -rn | tail -1 `";
	f << std::endl;
	f << "echo \"CPUs: $cpus\"" << std::endl;

	f << "pandda.analyse data_dirs='*' "
	<< "pdb_style='final.pdb' mtz_style='final.mtz' "
	<< "cpus=$cpus "
	<< "high_res_increment=0.3 min_build_datasets=20 "
	<< "checks.all_data_are_valid_values=None";
	
	f << std::endl;


	f.close();
}
