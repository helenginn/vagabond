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
#include <QToolBar>
#include <QApplication>
#include <QTabWidget>
#include <QActionGroup>
#include <QPushButton>
#include <QMessageBox>
#include <QMenu>
#include <QVariant>
#include <QVBoxLayout>
#include <QColor>
#include <QLabel>
#include <QPlainTextEdit>
#include <QFileDialog>
#include <QIcon>
#include <h3dsrc/Dialogue.h>
#include <iostream>
#include <iomanip>

#include "FileReader.h"
#include "Screen.h"
#include "SelectionWindow.h"
#include "CAlphaView.h"
#include "DisplaySettings.h"
#include "ClusterList.h"
#include "CorrelLabel.h"
#include "QuickAtoms.h"
#include "MtzFile.h"
#include "KeeperGL.h"
#include "GLPoint.h"
#include "Group.h"
#include "SelectMatrix.h"
#include "AxisScroll.h"
#include "MatrixView.h"
#include "HKLView.h"
#include "AveDiffraction.h"
#include "AveCSV.h"
#include "MtzFFTPtr.h"
#include "ColumnView.h"
#include "AveVectors.h"

#define TOOL_BAR_HEIGHT 50
#define TREE_VIEW_WIDTH 300
#define TAB_VIEW_WIDTH 800
#define RIGHT_VIEW_WIDTH 200

Screen::Screen(QWidget *widget) : QMainWindow(widget)
{
	_returnJourney = NULL;
	setGeometry(0, 0, 900, 600);

	_scale = -1;
	_storeHKL = make_mat4x4();
	_tabs = NULL;
	_correlLabel = NULL;
	_correlImage = NULL;
	_svdView = NULL;
	_ucView = NULL;
	_rView = NULL;
	_group = NULL;
	_newSel = NULL;
	_toggleDead = NULL;
	_markSele = NULL;
	_invertSele = NULL;
	_hkl = NULL;
	_hklKeeper = NULL;
	_cAlphaKeeper = NULL;
	_currIndex = 0;

	QWidget *window = new QWidget();
	QHBoxLayout *hout = new QHBoxLayout();
	window->setLayout(hout);
	QVBoxLayout *vout = new QVBoxLayout();

	_inputTree = new QTreeWidget(NULL);
	_inputTree->setMinimumSize(250, 0);
	_inputTree->setMaximumSize(300, 10000);
	vout->addWidget(_inputTree);

	_list = new ClusterList(_inputTree);
	_list->setScreen(this);

	_toolBar = new QToolBar(NULL);
	_toolBar->setMinimumSize(250, 0);
	_toolBar->setMaximumSize(300, 10000);
	vout->addWidget(_toolBar);
	addToolBar();

	hout->addLayout(vout);

	_tabs = new QTabWidget(NULL);
	_tabs->setMinimumSize(500, 500);
	hout->addWidget(_tabs);
	
	_side = new QWidget(NULL);
	_side->setMinimumSize(200, 0);
	_side->setMaximumSize(300, 10000);
	{
		QVBoxLayout *vout = new QVBoxLayout();
		_side->setLayout(vout);
	}
	hout->addWidget(_side);
	
	addSideButtons();
	
	connect(_list, &ClusterList::updateSelections,
	        this, &Screen::refreshSelection);

	setCentralWidget(window);
}

void Screen::updateToolbar(Group *grp)
{
	uses.amp->setChecked(false);
	uses.ca->setChecked(false);
	uses.uc->setChecked(false);
	uses.csv->setChecked(false);
	uses.vec->setChecked(false);

	switch (grp->getType())
	{
		case AveDiff:
		uses.amp->setChecked(true);
		break;
		
		case AveCA:
		uses.ca->setChecked(true);
		break;
		
		case AveUC:
		uses.uc->setChecked(true);
		break;
		
		case AveVec:
		uses.vec->setChecked(true);
		break;
		
		case AveComma:
		uses.csv->setChecked(true);
		break;
		
		default:
		break;
	}

	uses.top->setChecked(false);
	uses.orig->setChecked(false);
	uses.self->setChecked(false);

	switch (grp->getWhichGroup())
	{
		case GroupTop:
		uses.top->setChecked(true);
		break;
		
		case GroupOriginal:
		uses.orig->setChecked(true);
		break;
		
		case GroupMe:
		uses.self->setChecked(true);
		
		default:
		break;
	}
}

void Screen::addCSVSwitcher()
{
	if (AveCSV::csvCount() <= 1)
	{
		return;
	}

	_menu->addSeparator();
	
	QMenu *choice = _menu->findChild<QMenu *>("which_csv");
	
	if (choice == NULL)
	{
		choice = _menu->addMenu("Which CSV?");
		choice->setObjectName("which_csv");
	}

	choice->clear();

	for (size_t i = 0; i < AveCSV::csvCount(); i++)
	{
		std::string file = AveCSV::csvName(i);
		QAction *csv = choice->addAction(QString::fromStdString(file));
		connect(csv, &QAction::triggered, _list, 
		        [=]() { _list->switchCSV(i); });
	}
}

void Screen::addToolBar()
{
	QPushButton *a = new QPushButton("Set average");
	
	QMenu *m = new QMenu(a);
	_menu = m;
	uses.amp = m->addAction("Use amplitudes");
	uses.amp->setCheckable(true);
	connect(uses.amp, &QAction::triggered, _list, &ClusterList::recipAverage);

	uses.ca = m->addAction("Use C-alphas");
	uses.ca->setCheckable(true);
	connect(uses.ca, &QAction::triggered, _list, &ClusterList::pdbAverage);

	uses.uc = m->addAction("Use unit cell");
	uses.uc->setCheckable(true);
	connect(uses.uc, &QAction::triggered, 
	        _list, &ClusterList::unitCellAverage);
	
	uses.csv = m->addAction("Use pair-wise comparisons");
	uses.csv->setCheckable(true);
	connect(uses.csv, &QAction::triggered, 
	        _list, &ClusterList::csvAverage);
	
	uses.vec = m->addAction("Use vectors");
	uses.vec->setCheckable(true);
	connect(uses.vec, &QAction::triggered, 
	        _list, &ClusterList::vecAverage);

	m->addSeparator();

	uses.top = m->addAction("Use top level");
	uses.top->setCheckable(true);
	connect(uses.top, &QAction::triggered, _list, &ClusterList::topAverage);

	uses.orig = m->addAction("Use original");
	uses.orig->setCheckable(true);
	connect(uses.orig, &QAction::triggered, _list, 
	        &ClusterList::originalAverage);

	uses.self = m->addAction("Use this group");
	uses.self->setCheckable(true);
	connect(uses.self, &QAction::triggered, _list, &ClusterList::myAverage);

	a->setMenu(m);
	_toolBar->addWidget(a);

	_cluster = _toolBar->addAction("Cluster");
	connect(_cluster, &QAction::triggered, this, &Screen::clusterGroup);

	QIcon del = qApp->style()->standardIcon(QStyle::SP_TrashIcon);
	QAction *d = _toolBar->addAction(del, "Bin");
	connect(d, &QAction::triggered, this, &Screen::removeCluster);
}

void Screen::removeCluster()
{
	QTreeWidgetItem *item = _inputTree->currentItem();

	if (!Group::isGroup(item))
	{
		QMessageBox msgBox;
		msgBox.setText("You can only delete clusters.");
		msgBox.exec();
		return;
	}
	
	_list->removeCluster(static_cast<Group *>(item));
}

void Screen::clusterGroup()
{
	QTreeWidgetItem *item = _inputTree->currentItem();

	if (!Group::isGroup(item))
	{
		QMessageBox msgBox;
		msgBox.setText("Please choose the data set group for clustering.");
		msgBox.exec();
		return;
	}
	
	_list->cluster(static_cast<Group *>(item));
}

void Screen::binTab()
{
	for (size_t i = 0; i < _bin.size(); i++)
	{
		(*_bin[i])->hide();
		(*_bin[i])->deleteLater();
		*_bin[i] = 0;
	}

	_bin.clear();
}

void Screen::addCorrelImage(Group *ave)
{
	if (!ave->getCorrelMatrix())
	{
		return;
	}

	_correlImage = ave->getCorrelMatrix();
	_correlImage->updateSelection();

	CorrelLabel *l = new CorrelLabel(NULL, _correlImage, this);
	_correlLabel = l;
	_correlLabel->setScaledContents(true);
	_correlLabel->setFocusPolicy(Qt::StrongFocus);
	_tabs->addTab(l, "Correlation matrix");
	_bin.push_back((QWidget **)&_correlLabel);
}

void Screen::addPlotView(PlotView **view, Group *ave, 
                         std::string title, PlotType type)
{
	(*view) = new PlotView(type, this);
	_tabs->addTab((*view), QString::fromStdString(title));

	(*view)->setScreen(this);
	(*view)->setup(ave);

	_bin.push_back((QWidget **)&(*view));
}

void Screen::addColumnView(Group *ave)
{
	_columnView = new ColumnView(this, ave);
	_columnView->setList(_list);
	_tabs->addTab(_columnView, "Column explorer");
	_bin.push_back((QWidget **)&_columnView);

}

void Screen::addCAlphaView()
{
	_cAlpha = new QWidget(this);
	_tabs->addTab(_cAlpha, "C-alpha explorer");
	
	_cAlphaKeeper = new KeeperGL(_cAlpha);

	_bin.push_back((QWidget **)&_cAlphaKeeper);
	_bin.push_back((QWidget **)&_cAlpha);
}

void Screen::addHKLView(VagFFTPtr fft, std::string filename)
{
	if (!fft || fft->nn() == 0)
	{
		return;
	}

	_hkl = new QWidget(this);
	_tabs->addTab(_hkl, "HKL explorer");
	
	if (_scale < 0 || _scale != _scale)
	{
		double average = fft->averageAll();
		_scale = average * 5;
	}

	_hklKeeper = new KeeperGL(_hkl);
	_hklKeeper->addAxes();
	_hklKeeper->addHKLView(fft, _scale);
	_hklKeeper->setModelMatrix(_storeHKL);
	
	std::string ucInfo = filename + ": ";
	ucInfo += AveDiffraction::unitCellDesc(fft);
	
	_ucLabel = new QPlainTextEdit(QString::fromStdString(ucInfo), _hkl);
	_ucLabel->setReadOnly(true);
	_ucLabel->setGeometry(0, 0, 600, 60);
	_ucLabel->show();

	_bin.push_back((QWidget **)&_hklKeeper);
	_bin.push_back((QWidget **)&_ucLabel);
	_bin.push_back((QWidget **)&_hkl);
}

void Screen::displaySingle(MtzFFTPtr fft)
{
	_group = NULL;

	if (_tabs != NULL)
	{
		disconnect(_tabs, &QTabWidget::currentChanged, 
		           this, &Screen::changeIndex);
	}

	binTab();


	addHKLView(fft, fft->getMtzFile()->getFilename());
	addCAlphaView();
	vec3 centre = _list->topCluster()->getCentre();
	_cAlphaKeeper->addCAlphaView(fft->getMtzFile(), centre);

	/*
	_tabs->show();

	int top = 10;

	_toggleDead = new QPushButton("Toggle dead", this);
	_toggleDead ->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, top,
	                          RIGHT_VIEW_WIDTH - 20, 40);

	connect(_toggleDead, &QPushButton::clicked,
	        _list, &ClusterList::toggleDead);
	_toggleDead->show();
	
	_bin.push_back((QWidget **)&_toggleDead);
	*/

	_tabs->setCurrentIndex(_currIndex);
	connect(_tabs, &QTabWidget::currentChanged, 
	        this, &Screen::changeIndex);
}

void Screen::addColour(QString colour, QString display, QMenu *m)
{
	QAction *col = m->addAction(display);
	col->setProperty("colour", QColor(colour));
	connect(col, &QAction::triggered, this, &Screen::changeColour);
}

void Screen::displayResults(Group *ave)
{
	_group = ave;
	
	vec3 hklc = empty_vec3();
	if (_hklKeeper != NULL)
	{
		_storeHKL = _hklKeeper->getModel();
		hklc = _hklKeeper->getCentre();
	}
	
	bool stored = false;
	mat4x4 storeCAlpha = make_mat4x4();
	vec3 c = empty_vec3();

	if (_cAlphaKeeper != NULL)
	{
		stored = true;
		storeCAlpha = _cAlphaKeeper->getModel();
		c = _cAlphaKeeper->getCentre();
	}

	if (_tabs != NULL)
	{
		disconnect(_tabs, &QTabWidget::currentChanged, 
		           this, &Screen::changeIndex);
	}

	updateToolbar(ave);
	binTab();

	if (!ave->getCorrelMatrix())
	{
		return;
	}

	addCorrelImage(ave);
	addPlotView(&_svdView, ave, "SVD explorer", PlotSVD);
	std::string average = "Average";
	addHKLView(ave->getAverageFFT(), average);
	
	if (_hklKeeper != NULL)
	{
		_hklKeeper->changeCentre(hklc);
	}

	if (QuickAtoms::loadedAtoms())
	{
		addCAlphaView();
		_cAlphaKeeper->addCAlphaView(ave);
		if (stored)
		{
			_cAlphaKeeper->setModelMatrix(storeCAlpha);
			_cAlphaKeeper->changeCentre(c);
		}
	}

	if (_hklKeeper)
	{
		addPlotView(&_ucView, ave, "Misc properties", PlotUnitCell);
	}
	
	if (AveVectors::hasVectorData())
	{
		addColumnView(ave);
	}

	_tabs->setCurrentIndex(_currIndex);
	_tabs->show();
	connect(_tabs, &QTabWidget::tabBarClicked, 
	        this, &Screen::refocus);
	connect(_tabs, &QTabWidget::currentChanged, 
	        this, &Screen::changeIndex);
	refocus(_currIndex);
	relinkPixmap();
}

void Screen::exportText()
{
	QMessageBox msgBox;
	msgBox.setText("How would you like the dataset names represented "\
	               "in the text file?");
	QAbstractButton *fullpath, *fullnames, *nicknames;
	fullpath = msgBox.addButton("Full file path", QMessageBox::ActionRole);
	fullnames = msgBox.addButton("Full file name", QMessageBox::ActionRole);
	nicknames = msgBox.addButton("Wildcard nickname", QMessageBox::ActionRole);
	msgBox.addButton("Cancel", QMessageBox::RejectRole);

	msgBox.exec();
	
	QAbstractButton *chosen = msgBox.clickedButton();
	
	if (chosen == fullpath)
	{
		_list->exportAll(ExportFullPath);
	}
	else if (chosen == fullnames)
	{
		_list->exportAll(ExportFilename);
	}
	else if (chosen == nicknames)
	{
		_list->exportAll(ExportNickname);
	}
}

void Screen::addSideButton(QWidget **buttPtr, std::string title, 
                           int *top)
{
	QLayout *l = _side->layout();
	(*buttPtr) = new QPushButton(QString::fromStdString(title), NULL);
	l->addWidget(*buttPtr);

	*top += 50;
}

void Screen::writeAverageMTZ()
{
	if (_group == NULL)
	{
		return;
	}

	QFileDialog *f = new QFileDialog(this, "Write MTZ", 
	                                 "Reflection file (*.mtz)");
	f->setAcceptMode(QFileDialog::AcceptSave);
	f->setFileMode(QFileDialog::AnyFile);
	f->setOptions(QFileDialog::DontUseNativeDialog);
	f->show();

    QStringList fileNames;

    if (f->exec())
    {
        fileNames = f->selectedFiles();
    }
    
    if (fileNames.size() < 1)
    {
		return;
    }

	_group->writeHKL(fileNames[0].toStdString());

}

void Screen::changeColour()
{
	if (!_group)
	{
		return;
	}

	QObject *obj = QObject::sender();
	QVariant v = obj->property("colour");
	QColor c = v.value<QColor>();

	double r = c.red() / 255.;
	double b = c.green() / 255.;
	double g = c.blue() / 255.;
	double a = c.alpha() / 255.;
	
	if (c == Qt::transparent)
	{
		a = -1;
	}

	_group->changeColour(r, g, b, a);
	refreshSelection();
}

void Screen::markSelection()
{
	if (_group == NULL)
	{
		return;
	}

	bool mark = !(_list->getLastAverage()->isMarked());
	_list->getLastAverage()->setMarked(mark);
	refreshSelection();
}

void Screen::killSelection()
{
	if (_group == NULL)
	{
		return;
	}

	bool dead = !(_list->getLastAverage()->isDead());
	_list->getLastAverage()->setDead(dead);
	refreshSelection();
}

void Screen::newSelection()
{
	if (_group == NULL)
	{
		return;
	}

	std::vector<MtzFFTPtr> mtzs = _group->getMtzsFromSelection();
	_list->makeGroup(mtzs);
}

void Screen::relinkPixmap()
{
	if (!_correlImage || !_correlLabel)
	{
		std::cout << "Something doesn't exist" << std::endl;
		return;
	}

	QImage i = _correlImage->scaled(_correlLabel->width(), 
	                                _correlLabel->height(),
	                                Qt::IgnoreAspectRatio);
	_correlLabel->setPixmap(QPixmap::fromImage(i));
}


void Screen::refreshSelection()
{
	if (_correlImage && _correlLabel)
	{
		_correlImage->updateSelection();
		relinkPixmap();
	}
	
	std::cout << "Refreshing." << std::endl;
	
	if (_svdView)
	{
		std::cout << "Repopulating the SVD" << std::endl;
		_svdView->keeper()->getPlot()->repopulate();
	}
	
	if (_hklKeeper)
	{
		_hklKeeper->getHKLView()->repopulate();
	}
	
	if (_cAlphaKeeper)
	{
		_cAlphaKeeper->getCAlphaView()->repopulate();
	}
	
	if (_rView)
	{
		_rView->keeper()->getPlot()->repopulate();
	}
	
	if (_ucView)
	{
		_ucView->keeper()->getPlot()->repopulate();
	}
	
	_list->updateColours();
	emit refreshed();
}

void Screen::changeIndex(int index)
{
	if (_currIndex >= 0)
	{
		_currIndex = index;
	}
}

void Screen::refocus(int index)
{
	if (index == 1 && _svdView != NULL)
	{
		_svdView->setFocus();
	}
	else if (index <= 0 && _correlLabel != NULL)
	{
		_correlLabel->setFocus();
	}
	else
	{
		setFocus();
	}
}


void Screen::saveImages()
{
	if (_cAlphaKeeper)
	{
		std::string fnAlpha = findNextFilename("cAlphas.png");
		_cAlphaKeeper->saveImage(fnAlpha);
	}

}

void Screen::collapsePositions()
{
	_group->collapseDatasets();
	refreshSelection();
}

void Screen::coverage()
{
	std::vector<MtzFFTPtr> results = _group->orderByCoverage();
	_list->makeGroup(results);
}

void Screen::reorder()
{
	_list->reorderMTZs();

}

void Screen::reindex()
{
	SelectMatrix *matrix = new SelectMatrix(NULL);
	matrix->setList(_list);
	matrix->show();
}

void Screen::reorderByFile()
{
	QFileDialog *f = new QFileDialog(this, "Choose reordering file", 
	                                 "Text file (*.txt)");
	f->setFileMode(QFileDialog::AnyFile);
	f->setOptions(QFileDialog::DontUseNativeDialog);
	f->show();

    QStringList fileNames;

    if (f->exec())
    {
        fileNames = f->selectedFiles();
    }
    
    if (fileNames.size() < 1)
    {
		return;
    }

	f->deleteLater();
	std::string filename = fileNames[0].toStdString();
	std::string contents = get_file_contents(filename);
	_list->loadClusters(contents);

}

void Screen::keyReleaseEvent(QKeyEvent *e)
{
	if (_cAlphaKeeper)
	{
		_cAlphaKeeper->keyReleaseEvent(e);
	}
}

void Screen::keyPressEvent(QKeyEvent *e)
{
	if (_cAlphaKeeper)
	{
		_cAlphaKeeper->keyPressEvent(e);
	}
}

void Screen::returnToSender()
{
	_returnJourney->finished();
}

void Screen::rotateDegrees()
{
	_svdView->rotate();
}

void Screen::displaySettings()
{
	DisplaySettings *ds = new DisplaySettings(NULL);
	ds->setList(_list);
	ds->show();

}

void Screen::plotSpin()
{
	std::string folder = openDialogue(this, "Image folder", "", true, true);

	if (folder == "")
	{
		return;
	}

	std::string pattern = folder + "/fr*.png";
	std::vector<std::string> files = glob(pattern);

	for (size_t i = 0; i < files.size(); i++)
	{
		remove(files[i].c_str());
	}

	FileReader::setOutputDirectory(folder);
	int count = 0;

	for (size_t i = 0; i < 360; i++)
	{
		std::string number = i_to_str(i);
		std::string zeros = std::string(5 - number.length(), '0');
		std::string filename = "fr_" + zeros + number + ".png";
		std::string path = FileReader::addOutputDirectory(filename);
		std::cout << path << std::endl;

		_svdView->rotate(-1);
		_svdView->keeper()->saveImage(path);
	}
}


void Screen::addSideButtons()
{
	int top = 10;

	addSideButton((QWidget **)&_newSel, "New group", &top);
	connect(_newSel, &QPushButton::clicked,
	        this, &Screen::newSelection);

	addSideButton((QWidget **)&_invertSele, "Invert selection", &top);
	connect(_invertSele, &QPushButton::clicked,
	        _list, &ClusterList::invertSelection);

	addSideButton((QWidget **)&_markSele, "Toggle marked cluster", &top);
	connect(_markSele, &QPushButton::clicked,
	        this, &Screen::markSelection);

	addSideButton((QWidget **)&_deadSele, "Toggle all dead", &top);
	connect(_deadSele, &QPushButton::clicked,
	        this, &Screen::killSelection);

	addSideButton((QWidget **)&_changeColour, "Recolour all", &top);
	{
		QMenu *m = new QMenu(_changeColour);

		QAction *col = m->addAction("None");
		col->setProperty("colour", QColor(Qt::transparent));

		addColour("aqua", "Aqua", m);
		addColour("skyblue", "Sky blue", m);
		addColour("cornflowerblue", "Cornflower blue", m);
		addColour("blueviolet", "Blue violet", m);
		addColour("mediumpurple", "Medium purple", m);
		addColour("purple", "Purple", m);
		addColour("forestgreen", "Forest green", m);
		addColour("orange", "Orange", m);
		addColour("darkorange", "Dark orange", m);
		addColour("grey", "Grey", m);
		addColour("orchid", "Orchid", m);
		addColour("darkred", "Dark red", m);
		addColour("springgreen", "Spring green", m);
		addColour("limegreen", "Lime green", m);

		_changeColour->setMenu(m);
	}

	/*
	addSideButton((QWidget **)&_changeData, "Change data", &top);
	QMenu *m = new QMenu(_changeData);

	QAction *act = m->addAction("Collapse molecules to selection");
	connect(act, &QAction::triggered, this, &Screen::collapsePositions);

	act = m->addAction("Reindex selection");
	connect(act, &QAction::triggered,this, &Screen::reindex);

	_changeData->setMenu(m);
	_changeData->show();

	addSideButton((QWidget **)&_coverage, "Coverage order", &top);
	connect(_coverage, &QPushButton::clicked,
	        this, &Screen::coverage);
	*/

	addSideButton((QWidget **)&_coverage, "Write average MTZ", &top);
	connect(_coverage, &QPushButton::clicked,
	        this, &Screen::writeAverageMTZ);

	addSideButton((QWidget **)&_reorder, "Reorder datasets", &top);
	QMenu *m = new QMenu(_reorder);

	QAction *act = m->addAction("... by marked");
	connect(act, &QAction::triggered, this, &Screen::reorder);
	act = m->addAction("... by file");
	connect(act, &QAction::triggered, this, &Screen::reorderByFile);

	_reorder->setMenu(m);
	_reorder->show();

	int bottom = height() - 50;

	QLayout *l = _side->layout();
	
	static_cast<QBoxLayout *>(l)->addStretch(1);

	_export = new QPushButton("Export", NULL);

	m = new QMenu(_export);
	QAction *a1 = m->addAction("Text files only");
	connect(a1, &QAction::triggered, this, &Screen::exportText);
	QAction *a2 = m->addAction("Prepare directories");
	connect(a2, &QAction::triggered, _list, &ClusterList::prepDirs);
	if (_returnJourney != NULL)
	{
		QAction *a1 = m->addAction("Return to sender");
		connect(a1, &QAction::triggered, this, &Screen::returnToSender);
	}

	_export->setMenu(m);
	
	bottom -= 50;

	_images = new QPushButton("Image view", NULL);

	m = new QMenu(_images);
	{
		QAction *a3 = m->addAction("Display settings");
		connect(a3, &QAction::triggered, this, &Screen::displaySettings);
		QAction *a2 = m->addAction("Rotate plot");
		connect(a2, &QAction::triggered, this, &Screen::rotateDegrees);
		QAction *a1 = m->addAction("Plot spin movie");
		connect(a1, &QAction::triggered, this, &Screen::plotSpin);
		QAction *a4 = m->addAction("Export C-alpha image");
		connect(a4, &QAction::triggered, this, &Screen::saveImages);
	}

	_images->setMenu(m);
	_images->show();

	l->addWidget(_images);
	l->addWidget(_export);

	bottom -= 50;

}
