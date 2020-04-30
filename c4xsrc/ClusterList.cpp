// Fuck COV
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

#include <QTreeWidget>
#include <QApplication>
#include <QMessageBox>
#include <QThread>
#include <QKeyEvent>
#include <QPushButton>
#include <iostream>
#include <fstream>

#include "FileReader.h"
#include "ClusterList.h"
#include <libsrc/PDBReader.h>
#include "MtzFile.h"
#include "MtzFFT.h"
#include "Averager.h"
#include "Screen.h"

#include <libsrc/FFT.h>

ClusterList::ClusterList(QTreeWidget *widget)
{
	_sqlInput = false;
	_selectMode = false;
	_removeMode = false;
	_lastAverage = NULL;
	_res = 3.5;
	_worker = NULL;
	_widget = widget;
	_widget->setHeaderLabel("Dataset groups");
	connect(_widget, &QTreeWidget::currentItemChanged,
	        this, &ClusterList::displayResults);
	connect(this, &ClusterList::updateSelections,
	        this, &ClusterList::updateColours);
}

ClusterList::~ClusterList()
{

}

void ClusterList::load(std::vector<DatasetPath> paths)
{
	_sqlInput = false;
	_paths = paths;
	loadFiles();
}

void ClusterList::setFiles(std::vector<std::string> files)
{
	for (size_t i = 0; i < files.size(); i++)
	{
		DatasetPath path;
		path.refinement_id = -1;
		path.mtz_path = files[i];
		std::string base = getBaseFilenameWithPath(files[i]);
		std::string pdb = base + ".pdb";
		path.pdb_path = pdb;

		_paths.push_back(path);
	}
}

bool ClusterList::loadFiles()
{
	if (_sqlInput)
	{
		getFromDatabase();
		return true;
	}

	Averager *ave = new Averager(_widget);
	_clusters.push_back(ave);
	ave->setMaxResolution(_res);

	for (size_t i = 0; i < _paths.size(); i++)
	{
		MtzFile *file = new MtzFile(_paths[i].mtz_path);
		_files.push_back(file);

		DiffractionMtzPtr mtz = DiffractionMtzPtr(new DiffractionMtz());
		mtz->setNeedsRfree(false);
		mtz->setResLimit(_res);
		mtz->setFilename(_paths[i].mtz_path);
		file->setRefinementID(_paths[i].refinement_id);

		try
		{
			mtz->load();
		}
		catch (Shouter *s)
		{
			std::cout << "Ignoring mtz " << _paths[i].mtz_path << 
			" due to problems" << std::endl;
			continue;
		}

		CrystalPtr crystal;

		std::string pdb = _paths[i].pdb_path;

		if (file_exists(pdb))
		{
			PDBReader reader;
			reader.setFilename(pdb);
			reader.ignoreAtomsExcept("CA");
			crystal = reader.getCrystal();
			
			file->setRWork(reader.getRWork());
			file->setRFree(reader.getRFree());
		}

		ave->addMtz(mtz, file, crystal);
	}

	if (ave->mtzCount() == 0)
	{
		QMessageBox msgBox;
		msgBox.setText(tr("No data sets have been successfully loaded.\n"
		                  "Please pass MTZ files in on the command"
		                  " line."));
		QPushButton quit("Quit");
		msgBox.setDefaultButton(&quit);
		msgBox.exec();
		return false;
	}

	_widget->setCurrentItem(ave);
	return true;
}

void ClusterList::removeCluster(Averager *ave)
{
	if (ave->isMarked())
	{
		QMessageBox msgBox;
		msgBox.setText(tr("Please unmark this cluster before removing it."));
		msgBox.exec();
		return;
	}
	
	std::cout << "Removing cluster " <<
	ave->text(0).toStdString() << std::endl;
	int index = _widget->indexOfTopLevelItem(ave);
	
	if (index == 0)
	{
		QMessageBox msgBox;
		msgBox.setText(tr("Will not remove top level cluster."));
		msgBox.exec();
		return;
	}

	_widget->takeTopLevelItem(index);
}

void ClusterList::average(Averager *item)
{
	if (_worker && _worker->isRunning())
	{
		return;
	}
	
	if (!_worker)
	{
		_worker = new QThread();
	}

	item->moveToThread(_worker);
	connect(this, SIGNAL(average()),
	        item, SLOT(performAverage()));

	connect(item, SIGNAL(resultReady()), 
	        this, SLOT(handleResults()));
	connect(item, SIGNAL(failed()), 
	        this, SLOT(handleError()));
	_worker->start();

	emit average();
}

void ClusterList::cluster(Averager *item)
{
	if (_worker && _worker->isRunning())
	{
		return;
	}

	if (!_worker)
	{
		_worker = new QThread();
	}

	item->moveToThread(_worker);
	connect(this, SIGNAL(cluster()),
	        item, SLOT(performCluster()));

	connect(item, SIGNAL(resultReady()), 
	        this, SLOT(handleResults()));
	connect(item, SIGNAL(failed()), 
	        this, SLOT(handleError()));
	_worker->start();

	emit cluster();
}

Averager *ClusterList::topCluster()
{
	QTreeWidgetItem *item = _widget->topLevelItem(0);

	if (!Averager::isAverager(item))
	{
		return NULL;
	}

	Averager *obj = static_cast<Averager *>(item);
	return obj;
}

void ClusterList::exportAll()
{
	std::string fnall = findNextFilename("all_clusters.txt");
	std::string fnmarked = findNextFilename("marked_clusters.txt");
	std::string fndead = findNextFilename("dead_sets.txt");
	std::ofstream fall;
	std::ofstream fmarked;
	fall.open(fnall);
	fmarked.open(fnmarked);
	
	for (int i = 0; i < _widget->topLevelItemCount(); i++)
	{
		QTreeWidgetItem *item = _widget->topLevelItem(i);
		
		if (!Averager::isAverager(item))
		{
			continue;
		}
		
		Averager *obj = static_cast<Averager *>(item);
		obj->writeToStream(fall, false);
		obj->writeToStream(fmarked, true);
	}

	fall.close();
	fmarked.close();
	
	std::ofstream fdead;
	fdead.open(fndead);
	for (size_t i = 0; i < _files.size(); i++)
	{
		if (_files[i]->isDead())
		{
			fdead << _files[i]->getFilename() << std::endl;
		}
	}

	fdead.close();

	std::string m = "Written all clusters to " + fnall + "\n";
	m += "Written marked clusters only to " + fnmarked + "\n";
	m += "Written dead data sets to " + fndead + "\n";

	QMessageBox msgBox;
	msgBox.setText(QString::fromStdString(m));
	msgBox.exec();
}

void ClusterList::keyPressEvent(QKeyEvent *event)
{
	if (event->key() == Qt::Key_Shift)
	{
		_selectMode = true;
	}
	if (event->key() == Qt::Key_Control)
	{
		_removeMode = true;
	}
}

void ClusterList::keyReleaseEvent(QKeyEvent *event)
{
	if (event->key() == Qt::Key_Shift)
	{
		_selectMode = false;
	}
	if (event->key() == Qt::Key_Control)
	{
		_removeMode = false;
	}
}

void ClusterList::displayResults()
{
	QTreeWidgetItem *item = _widget->currentItem();

	if (!item)
	{
		return;
	}
	
	if (Averager::isAverager(item))
	{
		Averager *obj = static_cast<Averager *>(item);
		_lastAverage = obj;
		_screen->displayResults(obj);
	}
	else if (MtzFFT::isMtzFFT(item))
	{
		MtzFFT *obj = static_cast<MtzFFT *>(item);
		MtzFFTPtr vag = obj->shared_from_this();
		_screen->displaySingle(vag);
	}
}

void ClusterList::handleResults()
{
	std::cout << "Job complete." << std::endl;
	_worker = NULL;
	
	clearSelection();
	Averager *obj = static_cast<Averager *>(QObject::sender());
	_screen->displayResults(obj);
	disconnect(this, SIGNAL(average()), nullptr, nullptr);
	disconnect(this, SIGNAL(cluster()), nullptr, nullptr);
}

void ClusterList::handleError()
{
	Averager *obj = static_cast<Averager *>(QObject::sender());
	std::string e = obj->getError();

	QMessageBox msgBox;
	msgBox.setText(QString::fromStdString(e));
	msgBox.exec();
}

std::string isCommand(std::string command, std::string start)
{
	if (command.substr(0, start.length()) == start)
	{
		return command.substr(start.length(), std::string::npos);
	}
	
	return "";
}

void ClusterList::setCommands(std::vector<std::string> commands)
{
	_commands = commands;
	
	for (size_t i = 0; i < _commands.size(); i++)
	{
		std::string value = isCommand(_commands[i], "--max-res=");
		if (value.length())
		{
			_res = atof(value.c_str());
		}

		value = isCommand(_commands[i], "--sql=");
		to_lower(value);
		if (value == "yes" || value == "true")
		{
			std::cout << "Using SQL Input" << std::endl;
			_sqlInput = true;
		}
	}
}

void ClusterList::invertSelection()
{
	for (size_t i = 0; i < _lastAverage->mtzCount(); i++)
	{
		MtzFile *file = _lastAverage->getMtz(i)->getMtzFile();
		bool select = file->isSelected();
		file->setSelected(!select);
	}
	
	emit updateSelections();
}

void ClusterList::clearSelection()
{
	for (size_t i = 0; i < _files.size(); i++)
	{
		MtzFile *file = _files[i];
		file->setSelected(false);
	}

	emit updateSelections();
}

void ClusterList::updateColours()
{
	for (size_t i = 0; i < _clusters.size(); i++)
	{
		Averager *ave = _clusters[i];
		
		for (size_t j = 0; j < ave->mtzCount(); j++)
		{
			ave->getMtz(j)->updateText();
		}
	}

}

void ClusterList::topAverage()
{
	QTreeWidgetItem *item = _widget->topLevelItem(0);

	if (!item || !Averager::isAverager(item))
	{
		return;
	}

	Averager *ave = static_cast<Averager *>(item);
	_lastAverage->copyFromAverage(ave);
}

void ClusterList::resetAverage()
{
	_lastAverage->useOriginalAverage();
}

void ClusterList::pdbAverage()
{
	_lastAverage->setType(AveCAlpha);
	average(_lastAverage);
}

void ClusterList::makeGroup(std::vector<MtzFFTPtr> mtzs, bool withAve)
{
	if (mtzs.size() == 0)
	{
		QMessageBox msgBox;
		msgBox.setText(tr("No data sets have been selected."));
		msgBox.exec();
		return;
	}

	Averager *ave = new Averager(_widget);
	
	if (withAve && _lastAverage)
	{
		ave->copyFromAverage(_lastAverage);
	}
	
	_clusters.push_back(ave);
	ave->setMaxResolution(_res);

	for (size_t i = 0; i < mtzs.size(); i++)
	{
		MtzFFTPtr mtz = mtzs[i];
		ave->addMtz(mtz);
	}
	ave->updateText();
	
	_widget->addTopLevelItem(ave);
}

void ClusterList::toggleDead()
{
	QTreeWidgetItem *item = _widget->currentItem();

	if (!item || !MtzFFT::isMtzFFT(item))
	{
		return;
	}

	MtzFFT *obj = static_cast<MtzFFT *>(item);
	bool dead = obj->getMtzFile()->isDead();
	obj->getMtzFile()->setDead(!dead);

	updateSelections();
}
