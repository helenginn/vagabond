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

#include <libsrc/PDBReader.h>
#include <libsrc/Crystal.h>
#include "Group.h"
#include "FileReader.h"
#include "ClusterList.h"
#include "Output.h"
#include "AveCSV.h"
#include "MtzFile.h"
#include "MtzFFT.h"
#include "Screen.h"
#include "FolderInput.h"
#include "QuickAtoms.h"

#include <libsrc/FFT.h>

ClusterList::ClusterList(QTreeWidget *widget)
{
	_onlyLoad = false;
	_max = 0;
	_skip = 0;
	_streamInput = false;
	_sqlInput = false;
	_csvInput = false;
	_selectMode = false;
	_removeMode = false;
	_lastAverage = NULL;
	_res = 3.5;
	_worker = NULL;
	_widget = widget;
	_widget->setHeaderLabel("Dataset groups");
	connect(_widget, &QTreeWidget::itemClicked,
	        this, &ClusterList::displayResults);
	connect(_widget, &QTreeWidget::currentItemChanged,
	        this, &ClusterList::selectedResults);
	connect(this, &ClusterList::updateSelections,
	        this, &ClusterList::updateColours);

	_widget->setContextMenuPolicy(Qt::CustomContextMenu);
	connect(_widget, &QTreeWidget::customContextMenuRequested,
	        this, &ClusterList::prepareMenu);
}

void ClusterList::prepareMenu(const QPoint &p)
{
	QTreeWidgetItem *item = _widget->itemAt(p);

	if (item && !Group::isGroup(item))
	{
		MtzFFT *fft = static_cast<MtzFFT *>(item);
		MtzFile *file = fft->getMtzFile();
		file->flipSelected();
	}
	
	updateSelections();
}

ClusterList::~ClusterList()
{

}

void ClusterList::load(std::vector<DatasetPath> paths)
{
	_sqlInput = false;
	_csvInput = false;
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
		path.pandda_mtz = files[i];
		path.metadata = getBaseFilename(files[i]);
		std::string base = getBaseFilenameWithPath(files[i]);
		std::string pdb = base + ".pdb";
		path.pdb_path = pdb;

		_paths.push_back(path);
	}
}

void ClusterList::getFromFolders()
{
	FolderInput *input = new FolderInput(NULL);
	input->setList(this);
}

void ClusterList::getFromCSV(AveCSV *csv)
{
	csv->setList(this);
	csv->load();
}

void ClusterList::addCSVSwitcher()
{
	_screen->addCSVSwitcher();
}

void ClusterList::getFromCSV(std::string csv)
{
	std::vector<std::string> csvs = split(csv, ',');
	AveCSV *input = new AveCSV(NULL);
	input->setList(this);

	for (size_t i = 0; i < csvs.size(); i++)
	{
		input->setFilename(csvs[i]);
		input->load();
	}

	input->preparePaths();
	_csvGroup = input;
	_screen->addCSVSwitcher();
}

void ClusterList::cycleCSV(bool forward)
{
	int current = AveCSV::currentChoice();
	current = (current + (forward ? 1 : -1)) % AveCSV::csvCount();
	
	AveCSV::setChosen(current);
	csvAverage();
	_screen->clusterGroup();
}

void ClusterList::switchCSV(int c)
{
	AveCSV::setChosen(c);
	csvAverage();
}

void ClusterList::loadFromMultistate(std::string pdb)
{
	if (!file_exists(pdb))
	{
		std::cout << "File does not exist: " << pdb << std::endl;
	}

	Group *grp = new Group(NULL);
	std::string contents = get_file_contents(pdb);
	int state;
	
	{
		std::string name = pdb + "_state_" + i_to_str(state);
		MtzFile *file = new MtzFile(name);
		file->setPanddaName(name);
		file->setPdbPath(name);
		file->setMetadata(name);
		file->setRefinementID(name);
		_files.push_back(file);

		DiffractionMtzPtr mtz = DiffractionMtzPtr(new DiffractionMtz());
		mtz->setFilename(name);

		CrystalPtr crystal;

		std::string contents = "";
		PDBReader reader;
		reader.setContents(contents);
		reader.ignoreAtomsExcept("CA");
		crystal = reader.getCrystal();
			
		file->setCrystal(crystal);
		grp->addMtz(mtz, file);
	}

}

bool ClusterList::loadFiles()
{
	if (_sqlInput)
	{
		getFromDatabase();
		return true;
	}
	if (_streamInput)
	{
		getFromStream();
		return true;
	}
	else if (_csvInput)
	{
		getFromCSV(_csv);
		return true;
	}

	Group *grp = new Group(NULL);
	grp->setMaxResolution(_res);

	for (size_t i = 0; i < _paths.size(); i++)
	{
		MtzFile *file = new MtzFile(_paths[i].mtz_path);
		_files.push_back(file);

		DiffractionMtzPtr mtz = DiffractionMtzPtr(new DiffractionMtz());
		mtz->setNeedsRfree(false);
		mtz->setResLimit(_res);
		mtz->setFilename(_paths[i].mtz_path);

		file->setPanddaName(_paths[i].pandda_mtz);
		file->setPdbPath(_paths[i].pdb_path);
		file->setMetadata(_paths[i].metadata);
		file->setRefinementID(_paths[i].refinement_id);
		file->setCifPath(_paths[i].cif_path);
		file->setLigPath(_paths[i].lig_path);

		try
		{
			if (_paths[i].mtz_path.length())
			{
				mtz->load();
			}
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
			
			file->setCrystal(crystal);
			file->setRWork(reader.getRWork());
			file->setRFree(reader.getRFree());
		}

		grp->addMtz(mtz, file);
	}

	if (grp->mtzCount() == 0)
	{
		getFromFolders();
		return true;
	}

	grp->setTopGroup();
	_widget->addTopLevelItem(grp);
	_clusters.push_back(grp);
	_widget->setCurrentItem(grp);
	grp->performAverage();
	
	if (_preload.length())
	{
		std::string contents = get_file_contents(_preload);
		loadClusters(contents);
	}

	return true;
}

void ClusterList::removeCluster(Group *ave)
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

void ClusterList::average(Group *item)
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

void ClusterList::cluster(Group *item)
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

Group *ClusterList::topCluster()
{
	QTreeWidgetItem *item = _widget->topLevelItem(0);

	if (!Group::isGroup(item))
	{
		return NULL;
	}

	Group *obj = static_cast<Group *>(item);
	return obj;
}

void ClusterList::exportAll(ExportType type)
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
		
		if (!Group::isGroup(item))
		{
			continue;
		}
		
		Group *obj = static_cast<Group *>(item);
		obj->writeToStream(fall, type, false);
		obj->writeToStream(fmarked, type, true);
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
	m += "\nPlease cite cluster4x if you use it during your research!\n";
	m += "\nGinn, H. M. (2020). Pre-clustering data sets using cluster4x improves the signal-to-noise ratio of high-throughput crystallography drug-screening analysis. Acta Crystallographica Section D: Structural Biology, 76(11).\n";

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

void ClusterList::selectedResults()
{
	QTreeWidgetItem *item = _widget->currentItem();

	if (!item)
	{
		return;
	}
	
	if (Group::isGroup(item))
	{
		Group *obj = static_cast<Group *>(item);
		_lastAverage = obj;
	}
}

void ClusterList::displayResults()
{
	QTreeWidgetItem *item = _widget->currentItem();

	if (!item)
	{
		return;
	}
	
	if (Group::isGroup(item))
	{
		Group *obj = static_cast<Group *>(item);
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
	
	Group *obj = static_cast<Group *>(QObject::sender());
	obj->moveToThread(QThread::currentThread());

	_worker = NULL;

	_screen->displayResults(obj);
	clearSelection();

	disconnect(this, SIGNAL(average()), nullptr, nullptr);
	disconnect(this, SIGNAL(cluster()), nullptr, nullptr);

	disconnect(obj, SIGNAL(resultReady()), this, SLOT(handleResults()));
	disconnect(obj, SIGNAL(failed()), this, SLOT(handleError()));
}

void ClusterList::handleError()
{
	Group *obj = static_cast<Group *>(QObject::sender());
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

		value = isCommand(_commands[i], "--skip=");
		if (value.length())
		{
			_skip = atoi(value.c_str());
		}

		value = isCommand(_commands[i], "--max-sets=");
		if (value.length())
		{
			_max = atoi(value.c_str());
		}

		value = isCommand(_commands[i], "--geom=");
		if (value.length() > 0)
		{
			_geom = value;
		}

		value = isCommand(_commands[i], "--spg=");
		if (value.length() > 0)
		{
			_spg = value;
		}

		value = isCommand(_commands[i], "--stream=");
		if (value.length() > 0)
		{
			_streamInput = true;
			_stream = value;
		}

		value = isCommand(_commands[i], "--sql=");
		to_lower(value);
		if (value == "yes" || value == "true")
		{
			std::cout << "Using SQL Input" << std::endl;
			_sqlInput = true;
		}

		value = isCommand(_commands[i], "--csv=");
		if (value.length() > 0)
		{
			std::cout << "Using CSV Input" << std::endl;
			_csvInput = true;
			_csv = value;
		}

		value = isCommand(_commands[i], "--preload=");
		if (value.length() > 0)
		{
			std::cout << "Using preloaded groups" << std::endl;
			_preload = value;
		}

		value = isCommand(_commands[i], "--target-id=");
		if (value.length() > 0)
		{
			if (value.find("'") != std::string::npos)
			{
				std::cout << "Remove that quotation mark "
				"from your target-id - danger - not using" << std::endl;
			}
			else
			{
				std::cout << "Setting cluster target ID for SQL" << std::endl;
				_targetID = value;
			}
		}


		value = isCommand(_commands[i], "--only-load=");
		if (value.length() > 0)
		{
			std::cout << "Using preloaded groups" << std::endl;
			_preload = value;
			_onlyLoad = true;
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
		Group *ave = _clusters[i];
		
		for (size_t j = 0; j < ave->mtzCount(); j++)
		{
			ave->getMtz(j)->updateText();
		}
	}

}

void ClusterList::topAverage()
{
	_lastAverage->useAverageGroup(GroupTop);
	_screen->updateToolbar(_lastAverage);
}

void ClusterList::originalAverage()
{
	_lastAverage->useAverageGroup(GroupOriginal);
	_screen->updateToolbar(_lastAverage);
}

void ClusterList::myAverage()
{
	_lastAverage->useAverageGroup(GroupMe);
	_screen->updateToolbar(_lastAverage);
}

void ClusterList::pdbAverage()
{
	_lastAverage->useAverageType(AveCA);
	_screen->updateToolbar(_lastAverage);
}

void ClusterList::recipAverage()
{
	_lastAverage->useAverageType(AveDiff);
	_screen->updateToolbar(_lastAverage);
}

void ClusterList::csvAverage()
{
	_lastAverage->useAverageType(AveComma);
	_screen->updateToolbar(_lastAverage);
}

void ClusterList::unitCellAverage()
{
	_lastAverage->useAverageType(AveUC);
	_screen->updateToolbar(_lastAverage);
}

Group *ClusterList::makeGroup(std::vector<MtzFFTPtr> mtzs)
{
	if (mtzs.size() == 0)
	{
		QMessageBox msgBox;
		msgBox.setText(tr("No data sets have been selected."));
		msgBox.exec();
		return NULL;
	}

	Group *ave = new Group(_widget);
	
	if (_lastAverage)
	{
		ave->copyFromOriginal(_lastAverage);
		ave->useAverageGroup(GroupOriginal);
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
	return ave;
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

void ClusterList::prepDirs()
{
	Output o;

	int count = 0;
	for (int i = 0; i < _widget->topLevelItemCount(); i++)
	{
		QTreeWidgetItem *item = _widget->topLevelItem(i);
		if (!Group::isGroup(item))
		{
			continue;
		}

		Group *obj = static_cast<Group *>(item);
		count += o.prepCluster(obj);
	}

	QMessageBox msgBox;
	std::string str = ("Set up directories for " + 
	                   i_to_str(count) + " clusters.");
	str += "\nPlease cite cluster4x if you use it during your research!\n";
	str += "\nGinn, Acta Crystallographica D. 2020 (accepted).\n";
	str += "(full citation TBD...)\n";

	msgBox.setText(QString::fromStdString(str));
	msgBox.exec();
}

void ClusterList::reorderMTZs()
{
	std::vector<MtzFFTPtr> newList;
	std::vector<MtzFFTPtr> mtzs = Group::topGroup()->mtzs();

	for (int i = 0; i < _widget->topLevelItemCount(); i++)
	{
		QTreeWidgetItem *item = _widget->topLevelItem(i);
		if (!Group::isGroup(item))
		{
			continue;
		}

		Group *obj = static_cast<Group *>(item);
		if (!obj->isMarked())
		{
			continue;
		}
		
		std::vector<MtzFFTPtr> marked = obj->mtzs();
		
		for (size_t j = 0; j < marked.size(); j++)
		{
			newList.push_back(marked[j]);
		}
	}
	
	makeGroup(newList);
}

void ClusterList::setReindexMatrix(mat3x3 reindex, vec3 trans)
{
	mat3x3 inv = mat3x3_inverse(reindex);
	std::cout << mat3x3_desc(inv) << std::endl;

	QTreeWidgetItem *item = _widget->currentItem();

	if (!item || !Group::isGroup(item))
	{
		return;
	}

	Group *grp = static_cast<Group *>(item);
	for (size_t i = 0; i < grp->mtzCount(); i++)
	{
		if (!grp->getMtz(i)->getMtzFile()->isSelected())
		{
			continue;
		}

		grp->getMtz(i)->reindex(reindex);

		std::cout << "Reindexing " << i << std::endl;
		grp->getMtzFile(i)->getCrystal()->reindexAtoms(inv, trans);
		grp->getMtzFile(i)->getQuickAtoms()->fetchAtoms();
	}

	std::cout << "Here..." << std::endl;
	_screen->refreshSelection();
}

void ClusterList::loadClusters(std::string contents)
{
	std::vector<std::string> lines = split(contents, '\n');
	Group *top = Group::topGroup();

	for (size_t i = 0; i < lines.size(); i++)
	{
		std::string line = lines[i];
		std::vector<std::string> labels = split(line, ',');
		
		std::vector<MtzFFTPtr> files;
		for (size_t j = 0; j < labels.size(); j++)
		{
			std::string label = labels[j];
			
			for (size_t k = 0; k < top->mtzCount(); k++)
			{
				if (top->getMtzFile(k)->metadata() == label)
				{
					files.push_back(top->getMtz(k));
					break;
				}
			}
		}
		
		makeGroup(files);
	}
}

