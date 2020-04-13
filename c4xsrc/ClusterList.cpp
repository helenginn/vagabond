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
#include <QMessageBox>
#include <QThread>
#include <iostream>
#include <fstream>

#include "FileReader.h"
#include "ClusterList.h"
#include "MtzFile.h"
#include "MtzFFT.h"
#include "Averager.h"
#include "Screen.h"

#include <libsrc/FFT.h>

ClusterList::ClusterList(QTreeWidget *widget)
{
	_lastAverage = NULL;
	_worker = NULL;
	_widget = widget;
	_widget->setHeaderLabel("Dataset groups");
	connect(_widget, &QTreeWidget::currentItemChanged,
	        this, &ClusterList::displayResults);
}

void ClusterList::loadFiles(std::vector<std::string> files)
{
	try
	{
		double res = 3.5;
		Averager *ave = new Averager(_widget);
		_clusters.push_back(ave);
		ave->setMaxResolution(res);

		for (size_t i = 1; i < files.size(); i++)
		{
			MtzFile *file = new MtzFile(files[i]);
			_files.push_back(file);

			DiffractionMtzPtr mtz = DiffractionMtzPtr(new DiffractionMtz());
			mtz->setResLimit(res);
			mtz->setFilename(files[i]);
			mtz->load();

			ave->addMtz(mtz, file);
		}
	}
	catch (Shouter *s)
	{
		s->shoutToStdOut();
	}
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

std::string ClusterList::findNextFilename(std::string file)
{
	std::string trial = file;
	int count = 0;

	while (true)
	{
		std::string test = i_to_str(count) + "_" + trial;

		if (!file_exists(test))
		{
			return test;
		}
		
		count++;
	}
}

void ClusterList::exportAll()
{
	std::string fnall = findNextFilename("all_clusters.txt");
	std::string fnmarked = findNextFilename("marked_clusters.txt");
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

	std::string m = "Written all clusters to " + fnall + "\n";
	m += "Written marked clusters only to " + fnmarked + "\n";

	QMessageBox msgBox;
	msgBox.setText(QString::fromStdString(m));
	msgBox.exec();
}

void ClusterList::displayResults()
{
	QTreeWidgetItem *item = _widget->currentItem();

	if (!item)
	{
		return;
	}
	
	if (!Averager::isAverager(item))
	{
		return;
	}

	Averager *obj = static_cast<Averager *>(item);
	
	_lastAverage = obj;
	_screen->displayResults(obj);
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

void ClusterList::setCommands(std::vector<std::string> commands)
{
	_commands = commands;
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

void ClusterList::makeGroup(std::vector<MtzFFTPtr> mtzs, bool withAve)
{
	double res = 3.5;
	Averager *ave = new Averager(_widget);
	
	if (withAve && _lastAverage)
	{
		VagFFTPtr fft = _lastAverage->getAverageFFT();
		ave->setAverageFFT(fft);
	}
	
	_clusters.push_back(ave);
	ave->setMaxResolution(res);

	std::cout << "We have " << mtzs.size() << std::endl;
	for (size_t i = 0; i < mtzs.size(); i++)
	{
		MtzFFTPtr mtz = mtzs[i];
		ave->addMtz(mtz);
	}
	ave->updateText();
	
	_widget->addTopLevelItem(ave);
}
