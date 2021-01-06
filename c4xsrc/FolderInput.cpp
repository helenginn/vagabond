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

#include "FolderInput.h"
#include "ClusterList.h"
#include <QLabel>
#include <stdexcept>
#include <QLineEdit>
#include <QPushButton>
#include <QMessageBox>
#include <hcsrc/FileReader.h>

std::string findMutable(std::string path, std::string glob)
{
	int endPath = path.length();
	for (int end = glob.length() - 1; end >= 0; end--)
	{
		if (glob[end] != '*' || glob[end] != '?')
		{
			endPath--;
		}
		else
		{
			break;
		}
	}
	endPath++;

	size_t start = glob.find("*");
	
	std::string mut = path.substr(start, endPath);
	return mut;
}

FolderInput::FolderInput(QWidget *widget) : QMainWindow(widget)
{
	setWindowTitle("Load datasets");
	setGeometry(300, 300, 500, 540);
	show();
	
	int top = 40;
	
	QLabel *l = new QLabel("Folder pattern", this);
	l->setGeometry(40, top, 150, 40);
	l->show();
	_bin.push_back(l);

	QLineEdit *line = new QLineEdit(this);
	line->setPlaceholderText("relative/path/to/dataset-x000");
	line->setGeometry(150, top, 260, 40);
	line->show();
	_folderLine = line;
	_bin.push_back(line);

	top += 40;
	
	l = new QLabel("MTZ style", this);
	l->setGeometry(40, top, 150, 40);
	l->show();
	_bin.push_back(l);

	line = new QLineEdit(this);
	line->setPlaceholderText("placeholder.mtz");
	line->setGeometry(150, top, 260, 40);
	line->show();
	_mtzLine = line;
	_bin.push_back(line);
	
	top += 40;
	
	l = new QLabel("PDB style", this);
	l->setGeometry(40, top, 150, 40);
	l->show();
	_bin.push_back(l);

	line = new QLineEdit(this);
	line->setPlaceholderText("placeholder.pdb");
	line->setGeometry(150, top, 260, 40);
	line->show();
	_pdbLine = line;
	_bin.push_back(line);
	
	top += 40;
	
	l = new QLabel("Ligand style", this);
	l->setGeometry(40, top, 150, 40);
	l->show();
	_bin.push_back(l);

	line = new QLineEdit(this);
	line->setPlaceholderText("ligand.cif (optional)");
	line->setGeometry(150, top, 260, 40);
	line->show();
	_cifLine = line;
	_bin.push_back(line);
	
	top += 40;

	line = new QLineEdit(this);
	line->setPlaceholderText("ligand.pdb (optional)");
	line->setGeometry(150, top, 260, 40);
	line->show();
	_ligLine = line;
	_bin.push_back(line);
	
	top += 80;
	
	l = new QLabel("Get clusters:", this);
	l->setGeometry(40, top, 150, 40);
	l->show();
	_bin.push_back(l);

	line = new QLineEdit(this);
	line->setPlaceholderText("0_all_clusters.txt (optional)");
	line->setGeometry(150, top, 260, 40);
	line->show();
	_premake = line;
	_bin.push_back(line);

	top += 80;
	
	
	QPushButton *p = new QPushButton("Load", this);
	p->setGeometry(330, top, 80, 40);
	p->show();
	_bin.push_back(p);
	
	connect(p, &QPushButton::clicked, this, &FolderInput::load);

	top += 80;
	setGeometry(300, 300, 500, top);
}

FolderInput::~FolderInput() 
{
	for (size_t i = 0; i < _bin.size(); i++)
	{
		_bin[i]->deleteLater();
	}
}

void FolderInput::load()
{
	std::string folderGlob = _folderLine->text().toStdString();
	std::string mtzGlob = _mtzLine->text().toStdString();
	std::string pdbGlob = _pdbLine->text().toStdString();
	std::string cifGlob = _cifLine->text().toStdString();
	std::string ligGlob = _ligLine->text().toStdString();
	
	size_t found = 0;
	size_t missing = 0;

	std::vector<std::string> results;
	try
	{
		results = glob(folderGlob);
	}
	catch (std::runtime_error &err)
	{
	}

	std::cout << "Found " << results.size() << " results." << std::endl;
	
	if (results.size() == 0)
	{
		QMessageBox msgBox;
		msgBox.setText(tr("Folder search path does not find any results.\n"
		                  "Please check your search path."));
		msgBox.exec();
		_folderLine->setFocus();
		return;
	}
	
	std::vector<DatasetPath> paths;
	
	for (size_t i = 0; i < results.size(); i++)
	{
		std::vector<std::string> findMtz;
		std::vector<std::string> findPdb;
		std::vector<std::string> findCif;
		std::vector<std::string> findLig;
		DIR *dir = opendir(results[i].c_str());

		if (dir)
		{
			closedir(dir);
		}
		else
		{
			std::cout << results[i] << " not a directory." << std::endl;
			continue;
		}

		try
		{
			std::string search = results[i] + "/" + mtzGlob;
			findMtz = glob(search);
		}
		catch (std::runtime_error &err)
		{
		}

		if (findMtz.size() > 1)
		{
			QMessageBox msgBox;
			msgBox.setText(QString::fromStdString("Search for MTZ yielded "
			               + i_to_str(findMtz.size()) + " results for " 
			               + results[i] + ".\n" +
			               "Please make your MTZ search path "
			               "more specific to only find one MTZ."));
			msgBox.exec();
			return;
		}

		try
		{
			std::string search = results[i] + "/" + pdbGlob;
			findPdb = glob(search);
		}
		catch (std::runtime_error &err)
		{
			std::cout << err.what();
		}

		try
		{
			std::string search = results[i] + "/" + ligGlob;
			findLig = glob(search);
		}
		catch (std::runtime_error &err)
		{
			std::cout << err.what();
		}

		try
		{
			std::string search = results[i] + "/" + cifGlob;
			findCif = glob(search);
		}
		catch (std::runtime_error &err)
		{
			std::cout << err.what();
		}


		if (findPdb.size() > 1)
		{
			QMessageBox msgBox;
			msgBox.setText(QString::fromStdString("Search for PDB yielded " 
			                  + i_to_str(findPdb.size()) + " results for " 
			                  + results[i] + ".\n"
			                  "Please make your PDB search path "
			                  "more specific to only find one PDB."));
			msgBox.exec();
			return;
		}
		
		if (findCif.size() > 1 || findLig.size() > 1)
		{
			std::cout << "Warning: CIF file has multiple matches" 
			<< " in " << results[i] << std::endl;
		}
		
		if (findPdb.size() > 0 && findMtz.size() > 0)
		{
			found++;
			DatasetPath path;
			path.mtz_path = findMtz[0];
			path.pandda_mtz = findMtz[0];
			path.pdb_path = findPdb[0];
			
			if (findCif.size() == 1)
			{
				path.cif_path = findCif[0];
			}
			
			if (findLig.size() == 1)
			{
				path.lig_path = findLig[0];
			}

			path.metadata = findMutable(results[i], folderGlob);
			path.refinement_id = path.metadata;
			paths.push_back(path);
			std::cout << "Found files for " << path.metadata << std::endl;
		}
		else
		{
			std::cout << results[i] << " missing files." << std::endl;
			missing++;
		}
	}
	
	
	QMessageBox msgBox;
	msgBox.setText(QString::fromStdString("Of " + i_to_str(results.size()) + " folders, "
	                  "found " + i_to_str(found) + " usable folders, and "
	                  + i_to_str(missing) + " with at least one missing "
	                  "MTZ or PDB file."));

	_list->load(paths);
	hide();
	
	if (_premake->text().length())
	{
		std::string file = _premake->text().toStdString();
		std::string contents = get_file_contents(file);
		_list->loadClusters(contents);
	}

	msgBox.exec();
	deleteLater();
}

