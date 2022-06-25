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
#include "MtzOptions.h"
#include <QLabel>
#include <QLineEdit>
#include <QCheckBox>
#include <QPushButton>
#include <stdexcept>
#include <QHBoxLayout>
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
	QWidget *window = new QWidget();
	QVBoxLayout *vout = new QVBoxLayout();
	window->setLayout(vout);
	window->setMinimumSize(400, 0);
	_sub = NULL;

	setWindowTitle("Load datasets");
	
	{
		QHBoxLayout *hout = new QHBoxLayout();
		QLabel *l = new QLabel("Folder pattern", NULL);
		hout->addWidget(l);

		QLineEdit *e = new QLineEdit(NULL);
		e->setPlaceholderText("relative/path/to/dataset-x000");
		_folderLine = e;
		hout->addWidget(e);
		
		vout->addLayout(hout);
	}
	
	{
		QHBoxLayout *hout = new QHBoxLayout();
		QLabel *l = new QLabel("MTZ style", NULL);
		hout->addWidget(l);

		QLineEdit *e = new QLineEdit(this);
		e->setPlaceholderText("placeholder.mtz");
		_mtzLine = e;
		hout->addWidget(e);
		
		QPushButton *p = new QPushButton("Options...", NULL);
		connect(p, &QPushButton::clicked, this, &FolderInput::mtzOpts);
		hout->addWidget(p);
		
		vout->addLayout(hout);
	}
	
	{
		QHBoxLayout *hout = new QHBoxLayout();
		QLabel *l = new QLabel("PDB style", NULL);
		hout->addWidget(l);

		QLineEdit *e = new QLineEdit(this);
		e->setPlaceholderText("placeholder.pdb");
		_pdbLine = e;
		hout->addWidget(e);

		QCheckBox *c = new QCheckBox("Optional", this);
		c->setObjectName("optional_pdb");
		hout->addWidget(c);
		
		vout->addLayout(hout);
	}
	
	{
		QHBoxLayout *hout = new QHBoxLayout();
		QLabel *l = new QLabel("Ligand style", NULL);
		hout->addWidget(l);

		QVBoxLayout *v2out = new QVBoxLayout();
		{
			QLineEdit *e = new QLineEdit(this);
			e->setPlaceholderText("ligand.cif (optional)");
			_cifLine = e;
			v2out->addWidget(e);
		}

		{
			QLineEdit *e = new QLineEdit(this);
			e->setPlaceholderText("ligand.pdb (optional)");
			_ligLine = e;
			v2out->addWidget(e);
		}
		
		hout->addLayout(v2out);
		vout->addLayout(hout);
	}

	{
		QHBoxLayout *hout = new QHBoxLayout();
		QLabel *l = new QLabel("Get clusters" , NULL);
		hout->addWidget(l);

		QLineEdit *e = new QLineEdit(this);
		e->setPlaceholderText("0_all_clusters.txt (optional)");
		_premake = e;
		hout->addWidget(e);
		
		vout->addLayout(hout);
	}
	
	{
		QPushButton *p = new QPushButton("Load", NULL);
		connect(p, &QPushButton::clicked, this, &FolderInput::load);
		vout->addWidget(p);
	}
	
	setCentralWidget(window);
	show();
}

FolderInput::~FolderInput() 
{
	if (_sub != NULL)
	{
		delete _sub;
		_sub = NULL;
	}

}

void FolderInput::mtzOpts()
{
	if (_sub != NULL)
	{
		delete _sub;
		_sub = NULL;
	}

	MtzOptions *sub = new MtzOptions(NULL);
	sub->setList(_list);
	_sub = sub;
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
		results = glob_pattern(folderGlob);
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
	bool optional = findChild<QCheckBox *>("optional_pdb")->isChecked();
	
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
			findMtz = glob_pattern(search);
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
			findPdb = glob_pattern(search);
		}
		catch (std::runtime_error &err)
		{
			std::cout << err.what();
		}

		try
		{
			std::string search = results[i] + "/" + ligGlob;
			findLig = glob_pattern(search);
		}
		catch (std::runtime_error &err)
		{
			std::cout << err.what();
		}

		try
		{
			std::string search = results[i] + "/" + cifGlob;
			findCif = glob_pattern(search);
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
		
		if ((findPdb.size() > 0 || optional) && findMtz.size() > 0)
		{
			found++;
			DatasetPath path;
			path.mtz_path = findMtz[0];
			path.pandda_mtz = findMtz[0];
			if (findPdb.size() > 0)
			{
				path.pdb_path = findPdb[0];
			}
			
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

