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

#include "Input.h"
#include "DatasetPath.h"
#include "FolderInput.h"
#include "ClusterList.h"
#include <QVBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <h3dsrc/Dialogue.h>
#include <hcsrc/FileReader.h>

Input::Input(QWidget *parent) : QMainWindow(parent)
{
	QWidget *window = new QWidget();
	QVBoxLayout *vout = new QVBoxLayout();
	window->setLayout(vout);
	window->setMinimumSize(400, 0);

	setWindowTitle("Data entry type");

	{
		QLabel *l = new QLabel("Load crystallographic datasets");
		vout->addWidget(l);
	}

	{
		QPushButton *p = new QPushButton("... from search pattern", NULL);
		connect(p, &QPushButton::clicked, this, &Input::loadCrystals);
		vout->addWidget(p);
	}

	{
		QPushButton *p = new QPushButton("... from list", NULL);
		connect(p, &QPushButton::clicked, this, &Input::loadCrystalList);
		vout->addWidget(p);
	}

	{
		QLabel *l = new QLabel("Load data from CSV");
		vout->addWidget(l);
	}

	{
		QPushButton *p = new QPushButton("Samples with arbitrary data", NULL);
		connect(p, &QPushButton::clicked, this, &Input::loadArbitrary);
		vout->addWidget(p);
	}

	{
		QPushButton *p = new QPushButton("List of pair-wise comparisons", NULL);
		connect(p, &QPushButton::clicked, this, &Input::loadCSV);
		vout->addWidget(p);
	}

	setCentralWidget(window);
	show();
}

void Input::loadCrystals()
{
	FolderInput *input = new FolderInput(NULL);
	input->setList(_list);
	this->hide();
}

void Input::loadCrystalList()
{
	std::string filename;
	filename = openDialogue(this, "Load CSV", 
	                        "Comma-separated values (*.csv)",
	                        false);

	try
	{
		std::string f = get_file_contents(filename);
		std::vector<std::string> lines = split(f, '\n');
		
		std::vector<DatasetPath> paths;
		for (size_t i = 0; i < lines.size(); i++)
		{
			std::vector<std::string> bits = split(lines[i], '\n');
			bits.resize(3);
			for (size_t j = 0; j < bits.size(); j++)
			{
				trim(bits[j]);
			}
			DatasetPath path;
			path.refinement_id = -1;
			path.mtz_path = bits[0];
			path.pandda_mtz = bits[0];
			path.metadata = getBaseFilename(bits[0]);
			std::string base = getBaseFilenameWithPath(bits[0]);
			std::string pdb = base + ".pdb";
			path.pdb_path = pdb;
			if (bits[1].length())
			{
				path.pdb_path = bits[1];
			}
			if (bits[2].length())
			{
				path.lig_path = bits[2];
			}
			
			paths.push_back(path);
		}
		
		_list->load(paths);
		hide();
	}
	catch (int e)
	{
		QMessageBox msg;
		msg.setText("Could not open file.");
		msg.exec();
		return;
	}
}

void Input::loadCSV()
{
	std::string filename;
	filename = openDialogue(this, "Load CSV", 
	                        "Comma-separated values (*.csv)",
	                        false);

	try
	{
		_list->getFromCSV(filename);
		hide();
	}
	catch (int e)
	{
		QMessageBox msg;
		msg.setText("Could not open file.");
		msg.exec();
		return;
	}
}

void Input::loadArbitrary()
{
	std::string filename;
	filename = openDialogue(this, "Load data", 
	                        "Comma-separated values (*.csv)",
	                        false);

	try
	{
		_list->loadFromVectorList(filename);
		hide();
	}
	catch (int e)
	{
		QMessageBox msg;
		msg.setText("Could not open file.");
		msg.exec();
		return;
	}
}
