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
#include <QColor>
#include <QLabel>
#include <QPlainTextEdit>
#include <QIcon>
#include <iostream>
#include <iomanip>

#include "FileReader.h"
#include "Screen.h"
#include "SelectionWindow.h"
#include "CAlphaView.h"
#include "ClusterList.h"
#include "CorrelLabel.h"
#include "MtzFile.h"
#include "KeeperGL.h"
#include "GLPoint.h"
#include "Group.h"
#include "AxisScroll.h"
#include "MatrixView.h"
#include "AveDiffraction.h"
#include "AveCSV.h"
#include "MtzFFTPtr.h"

#define TOOL_BAR_HEIGHT 50
#define TREE_VIEW_WIDTH 300
#define TAB_VIEW_WIDTH 800
#define RIGHT_VIEW_WIDTH 200

Screen::Screen(QWidget *widget) : QMainWindow(widget)
{
	setGeometry(0, 0, 1200, 800);

	_scale = -1;
	_storeHKL = make_mat4x4();
	_storeCAlpha = make_mat4x4();
	mat4x4_mult_scalar(&_storeHKL, 5);
	mat4x4_mult_scalar(&_storeCAlpha, 0.1);
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

	_inputTree = new QTreeWidget(this);
	_inputTree->show();

	_list = new ClusterList(_inputTree);
	_list->setScreen(this);
	addToolBar();
	
	connect(_list, &ClusterList::updateSelections,
	        this, &Screen::refreshSelection);

	resizeEvent(NULL);
}

void Screen::resizeEvent(QResizeEvent *e)
{
	int w = width();
	int h = height();
	int tool_h = TOOL_BAR_HEIGHT;

	int tree_w = std::max((int)(w * 0.2), TREE_VIEW_WIDTH);
	_inputTree->setGeometry(0,0, tree_w, h - tool_h);

	_toolBar->setGeometry(0, height() - TOOL_BAR_HEIGHT, 
	                      TREE_VIEW_WIDTH, TOOL_BAR_HEIGHT);
	
	if (_tabs)
	{
		_tabs->setGeometry(tree_w, 0, width() - tree_w 
		                   - RIGHT_VIEW_WIDTH, height());
		if (_correlLabel)
		{
			_correlLabel->setGeometry(0, tool_h, _tabs->width(), 
			                          _tabs->height() - tool_h);

			relinkPixmap();
		}

		if (_hklKeeper)
		{
			_hklKeeper->setGeometry(0, 0, _tabs->width(), 
			                        _tabs->height());
			
			if (_ucLabel)
			{
				_ucLabel->setGeometry(0, 0, _tabs->width(), 40);
			}
		}

		if (_cAlphaKeeper)
		{
			_cAlphaKeeper->setGeometry(0, 0, _tabs->width(), 
			                           _tabs->height());
		}

		if (_svdView)
		{
			_svdView->resizeEvent(e);
		}

		if (_ucView)
		{
			_ucView->resizeEvent(e);
		}

		if (_rView)
		{
			_ucView->resizeEvent(e);
		}

		if (_newSel)
		{
			_newSel->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, 
			                     _newSel->y(), RIGHT_VIEW_WIDTH - 20, 40);

			_markSele->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
			                       _markSele->y(), RIGHT_VIEW_WIDTH - 20,
			                       40);

			_invertSele->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
			                         _invertSele->y(), RIGHT_VIEW_WIDTH -
			                         20, 40);

			_deadSele->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, 
			                       _deadSele->y(), RIGHT_VIEW_WIDTH - 20,
			                       40);
			_collapse->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
			                       _collapse->y(), RIGHT_VIEW_WIDTH - 20,
			                       40);
			_coverage->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
			                       _coverage->y(), RIGHT_VIEW_WIDTH - 20,
			                       40);
			_reorder->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
			                       _reorder->y(), RIGHT_VIEW_WIDTH - 20,
			                       40);

			_changeColour->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
			                       _changeColour->y(), RIGHT_VIEW_WIDTH - 20,
			                           40);

			_export->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
			                     h - 50, RIGHT_VIEW_WIDTH - 20, 40);

			_images->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
			                     h - 100, RIGHT_VIEW_WIDTH - 20, 40);

		}

		if (_toggleDead)
		{
			_toggleDead->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, 
			                         _toggleDead->y(), RIGHT_VIEW_WIDTH - 20,
			                         60);
		}
	}
}

void Screen::updateToolbar(Group *grp)
{
	uses.amp->setChecked(false);
	uses.ca->setChecked(false);
	uses.uc->setChecked(false);
	uses.csv->setChecked(false);

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

void Screen::addToolBar()
{
	_toolBar = new QToolBar(this);
	_toolBar->show();
	
	QPushButton *a = new QPushButton("Set average");
	
	QMenu *m = new QMenu(a);
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
	
	uses.csv = m->addAction("Use CSV file");
	uses.csv->setCheckable(true);
	connect(uses.csv, &QAction::triggered, 
	        _list, &ClusterList::csvAverage);

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
	delete _tabs;
	_tabs = NULL;
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

void Screen::addCAlphaView()
{
	_cAlpha = new QWidget(this);
	_tabs->addTab(_cAlpha, "C-alpha explorer");
	
	_cAlphaKeeper = new KeeperGL(_cAlpha);
	_cAlphaKeeper->setStoreMatrix(&_storeCAlpha);

	_bin.push_back((QWidget **)&_cAlphaKeeper);
	_bin.push_back((QWidget **)&_cAlpha);
}

void Screen::addHKLView(VagFFTPtr fft, std::string filename)
{
	if (!fft)
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
	_hklKeeper->setStoreMatrix(&_storeHKL);
	_hklKeeper->addHKLView(fft, _scale);
	
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
	_tabs = new QTabWidget(this);
	resizeEvent(NULL);


	addHKLView(fft, fft->getMtzFile()->getFilename());
	addCAlphaView();
	vec3 centre = _list->topCluster()->getCentre();
	_cAlphaKeeper->addCAlphaView(fft->getMtzFile(), centre);

	_tabs->show();

	int top = 10;

	_toggleDead = new QPushButton("Toggle dead", this);
	_toggleDead ->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, top,
	                          RIGHT_VIEW_WIDTH - 20, 40);

	connect(_toggleDead, &QPushButton::clicked,
	        _list, &ClusterList::toggleDead);
	_toggleDead->show();
	
	_bin.push_back((QWidget **)&_toggleDead);

	_tabs->setCurrentIndex(_currIndex);
	connect(_tabs, &QTabWidget::currentChanged, 
	        this, &Screen::changeIndex);
	resizeEvent(NULL);
}

void Screen::displayResults(Group *ave)
{
	_group = ave;

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

	_tabs = new QTabWidget(this);

	resizeEvent(NULL);

	addCorrelImage(ave);
	addPlotView(&_svdView, ave, "SVD explorer", PlotSVD);
	std::string average = "Average";
	addHKLView(ave->getAverageFFT(), average);

	addCAlphaView();
	_cAlphaKeeper->addCAlphaView(ave);

	addPlotView(&_ucView, ave, "Misc properties", PlotUnitCell);
//	addPlotView(&_rView, ave, "R factors", PlotRFactors);

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
	_deadSele->setDisabled(_group->isTopGroup());
	connect(_deadSele, &QPushButton::clicked,
	        this, &Screen::killSelection);

	addSideButton((QWidget **)&_changeColour, "Recolour all", &top);
	{
		QMenu *m = new QMenu(_changeColour);

		QAction *col = m->addAction("None");
		col->setProperty("colour", QColor(Qt::transparent));

		col = m->addAction("Aqua");
		col->setProperty("colour", QColor("aqua"));
		connect(col, &QAction::triggered, this, &Screen::changeColour);

		col = m->addAction("Sky blue");
		col->setProperty("colour", QColor("skyblue"));
		connect(col, &QAction::triggered, this, &Screen::changeColour);

		col = m->addAction("Blue violet");
		col->setProperty("colour", QColor("blueviolet"));
		connect(col, &QAction::triggered, this, &Screen::changeColour);

		col = m->addAction("Medium purple");
		col->setProperty("colour", QColor("mediumpurple"));
		connect(col, &QAction::triggered, this, &Screen::changeColour);

		col = m->addAction("Forest green");
		col->setProperty("colour", QColor("forestgreen"));
		connect(col, &QAction::triggered, this, &Screen::changeColour);

		col = m->addAction("Orchid");
		col->setProperty("colour", QColor("orchid"));
		connect(col, &QAction::triggered, this, &Screen::changeColour);
		
		_changeColour->setMenu(m);
	}

	addSideButton((QWidget **)&_collapse, "Collapse positions", &top);
	connect(_collapse, &QPushButton::clicked,
	        this, &Screen::collapsePositions);

	addSideButton((QWidget **)&_coverage, "Coverage order", &top);
	connect(_coverage, &QPushButton::clicked,
	        this, &Screen::coverage);

	addSideButton((QWidget **)&_reorder, "Reorder by marked", &top);
	connect(_reorder, &QPushButton::clicked,
	        this, &Screen::reorder);

	int bottom = height() - 50;

	_export = new QPushButton("Export", this);
	_export->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, bottom,
	                     RIGHT_VIEW_WIDTH - 20, 40);

	QMenu *m = new QMenu(_export);
	QAction *a1 = m->addAction("Text files only");
	connect(a1, &QAction::triggered, _list, &ClusterList::exportAll);
	QAction *a2 = m->addAction("Prepare directories");
	connect(a2, &QAction::triggered, _list, &ClusterList::prepDirs);

	_export->setMenu(m);
	_export->show();
	
	bottom -= 50;

	_images = new QPushButton("Export images", this);
	_images->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, bottom,
	                     RIGHT_VIEW_WIDTH - 20, 40);
	_images->show();
	connect(_images, &QPushButton::clicked,
	        this, &Screen::saveImages);

	bottom -= 50;

	_bin.push_back((QWidget **)&_export);
	_bin.push_back((QWidget **)&_images);

	_tabs->setCurrentIndex(_currIndex);
	_tabs->show();
	connect(_tabs, &QTabWidget::tabBarClicked, 
	        this, &Screen::refocus);
	connect(_tabs, &QTabWidget::currentChanged, 
	        this, &Screen::changeIndex);
	refocus(_currIndex);
	relinkPixmap();
	resizeEvent(NULL);
}

void Screen::addSideButton(QWidget **buttPtr, std::string title, 
                           int *top)
{
	(*buttPtr) = new QPushButton(QString::fromStdString(title), this);
	(*buttPtr)->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, *top,
	               RIGHT_VIEW_WIDTH - 20, 40);
	(*buttPtr)->show();
	_bin.push_back(buttPtr);

	*top += 50;
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
	bool mark = !(_group->isMarked());
	_group->setMarked(mark);
	refreshSelection();
}

void Screen::killSelection()
{
	bool dead = !(_group->isDead());
	_group->setDead(dead);
	refreshSelection();
}

void Screen::newSelection()
{
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
		resizeEvent(NULL);
		relinkPixmap();
	}
	
	if (_svdView)
	{
		_svdView->keeper()->getPlot()->repopulate();
	}
	
	if (_ucView)
	{
		_ucView->keeper()->getPlot()->repopulate();
	}
	
	if (_rView)
	{
		_rView->keeper()->getPlot()->repopulate();
	}
	
	if (_cAlphaKeeper)
	{
		_cAlphaKeeper->getCAlphaView()->repopulate();
	}
	
	_list->updateColours();
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
		resizeEvent(NULL);
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
	Group *topGroup = _list->topCluster();
	_group->collapseDatasets(topGroup);
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
