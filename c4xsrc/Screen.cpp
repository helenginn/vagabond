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
#include <QLabel>
#include <QPlainTextEdit>
#include <QIcon>
#include <iostream>

#include "FileReader.h"
#include "Screen.h"
#include "SelectionWindow.h"
#include "CAlphaView.h"
#include "ClusterList.h"
#include "CorrelLabel.h"
#include "MtzFile.h"
#include "KeeperGL.h"
#include "GLPoint.h"
#include "Averager.h"
#include "AxisScroll.h"
#include "MatrixView.h"
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
	_graph = NULL;
	_keeper = NULL;
	_scroll = NULL;
	_newSel = NULL;
	_toggleDead = NULL;
	_markSele = NULL;
	_unmarkSele = NULL;
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

		if (_graph)
		{
			int axh = 60;
			
			if (_keeper)
			{
				_keeper->setGeometry(0, axh, _tabs->width(), 
				                     _tabs->height() - axh);
			}
			
			if (_selection)
			{
				_selection->setGeometry(-2, axh - 2, _tabs->width() + 2, 
				                        _tabs->height() - axh);
			}
			
			if (_scroll)
			{
				_scroll->setGeometry(0, 0, _tabs->width(), axh);
			}

			if (_newSel)
			{
				_newSel->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, 
				                     _newSel->y(), RIGHT_VIEW_WIDTH - 20, 40);

				_markSele->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
				                       _markSele->y(), RIGHT_VIEW_WIDTH - 20,
				                       40);
				_unmarkSele->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
										 _unmarkSele->y(), RIGHT_VIEW_WIDTH -
										                   20, 40);

				_invertSele->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
				                         _invertSele->y(), RIGHT_VIEW_WIDTH -
				                         20, 40);

				_deadSele->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, 
				                       _deadSele->y(), RIGHT_VIEW_WIDTH - 20,
				                       40);
				_undeadSele->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
				                         _undeadSele->y(), RIGHT_VIEW_WIDTH -
				                         20, 40);

				_export->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
				                     h - 50, RIGHT_VIEW_WIDTH - 20, 40);

				_images->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
				                     h - 100, RIGHT_VIEW_WIDTH - 20, 40);
			}
		}

		if (_toggleDead)
		{
			_toggleDead->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, 
			                         _toggleDead->y(), RIGHT_VIEW_WIDTH - 20,
			                         60);
		}
	}
}

void Screen::addToolBar()
{
	_toolBar = new QToolBar(this);
	_toolBar->show();
	
	QPushButton *a = new QPushButton("Set average");
	
	QMenu *m = new QMenu(a);
	QAction *a1 = m->addAction("Recalculate, diffraction");
	connect(a1, &QAction::triggered, this, &Screen::averageGroup);
	QAction *a4 = m->addAction("Recalculate, C-alphas");
	connect(a4, &QAction::triggered, _list, &ClusterList::pdbAverage);
	QAction *a2 = m->addAction("Use original");
	connect(a2, &QAction::triggered, _list, &ClusterList::resetAverage);
	QAction *a3 = m->addAction("Use top level");
	connect(a3, &QAction::triggered, _list, &ClusterList::topAverage);
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

	if (!Averager::isAverager(item))
	{
		QMessageBox msgBox;
		msgBox.setText("You can only delete clusters.");
		msgBox.exec();
		return;
	}
	
	_list->removeCluster(static_cast<Averager *>(item));
}

void Screen::averageGroup()
{
	QTreeWidgetItem *item = _inputTree->currentItem();

	if (!Averager::isAverager(item))
	{
		QMessageBox msgBox;
		msgBox.setText("Please choose the data set group for averaging.");
		msgBox.exec();
		return;
	}
	
	Averager *ave = static_cast<Averager*>(item);
	ave->setType(AveDiffraction);
	_list->average(ave);
}

void Screen::clusterGroup()
{
	QTreeWidgetItem *item = _inputTree->currentItem();

	if (!Averager::isAverager(item))
	{
		QMessageBox msgBox;
		msgBox.setText("Please choose the data set group for clustering.");
		msgBox.exec();
		return;
	}
	
	_list->cluster(static_cast<Averager *>(item));
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

void Screen::addCorrelImage(Averager *ave)
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

void Screen::addAxisExplorer(Averager *ave)
{
	if (!ave->getCorrelMatrix())
	{
		return;
	}

	_graph = new QWidget(this);
	_tabs->addTab(_graph, "Axis explorer");

	_keeper = new KeeperGL(_graph);
	_keeper->addAxes();
	_keeper->addSVDPoints(ave);

	GLPoint *points = _keeper->getPoints();

	connect(points, &GLPoint::updateSelection,
	        this, &Screen::refreshSelection);

	_selection = new SelectionWindow(_graph, _keeper);
	_selection->setPoints(points);
	_selection->show();
	_selection->setFocusPolicy(Qt::StrongFocus);

	_scroll = new AxisScroll(_graph);
	_scroll->setAverager(ave);
	_scroll->setPoints(points);
	_scroll->makeLayout();

	_bin.push_back((QWidget **)&_graph);
	_bin.push_back((QWidget **)&_scroll);
	_bin.push_back((QWidget **)&_selection);
	_bin.push_back((QWidget **)&_keeper);
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
	
	std::vector<double> uc = fft->getUnitCell();
	std::string ucInfo;
	ucInfo += filename + ": ";
	ucInfo += fft->getSpaceGroup()->symbol_Hall;
	ucInfo += "; a, b, c =  ";
	for (int i = 0; i < 3; i++)
	{
		ucInfo += f_to_str(uc[i], 2) + "   ";
	}
	ucInfo += " Å; a, b, g =  ";
	for (int i = 3; i < 6; i++)
	{
		ucInfo += f_to_str(uc[i], 2) + "   ";
	}
	ucInfo += "°.";
	
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

void Screen::displayResults(Averager *ave)
{
	if (_tabs != NULL)
	{
		disconnect(_tabs, &QTabWidget::currentChanged, 
		           this, &Screen::changeIndex);
	}

	binTab();
	if (!ave->getCorrelMatrix())
	{
		return;
	}
	_tabs = new QTabWidget(this);

	resizeEvent(NULL);

	addCorrelImage(ave);
	addAxisExplorer(ave);
	std::string average = "Average";
	addHKLView(ave->getAverageFFT(), average);

	addCAlphaView();
	vec3 centre = _list->topCluster()->getCentre();
	_cAlphaKeeper->addCAlphaView(ave, centre);

	int top = 10;

	_newSel = new QPushButton("New group", this);
	_newSel->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, top,
	                     RIGHT_VIEW_WIDTH - 20, 40);
	connect(_newSel, &QPushButton::clicked,
	        this, &Screen::newSelection);
	_newSel->show();

	top += 50;

	_invertSele = new QPushButton("Invert selection", this);
	_invertSele->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, top,
	                         RIGHT_VIEW_WIDTH - 20, 40);
	_invertSele->show();
	connect(_invertSele, &QPushButton::clicked,
	        _list, &ClusterList::invertSelection);

	top += 50;

	_markSele = new QPushButton("Mark all in cluster", this);
	_markSele->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, top,
	                       RIGHT_VIEW_WIDTH - 20, 40);
	_markSele->show();
	connect(_markSele, &QPushButton::clicked,
	        this, &Screen::markSelection);

	top += 50;

	_unmarkSele = new QPushButton("Unmark datasets", this);
	_unmarkSele->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, top,
	                         RIGHT_VIEW_WIDTH - 20, 40);
	_unmarkSele->show();
	connect(_unmarkSele, &QPushButton::clicked,
	        this, &Screen::markSelection);

	top += 50;

	_deadSele = new QPushButton("Mark all dead/rubbish", this);
	_deadSele->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, top,
	                       RIGHT_VIEW_WIDTH - 20, 40);
	_deadSele->show();
	connect(_deadSele, &QPushButton::clicked,
	        this, &Screen::killSelection);

	top += 50;

	_undeadSele = new QPushButton("Unmark dead/rubbish", this);
	_undeadSele->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, top,
	                         RIGHT_VIEW_WIDTH - 20, 40);
	_undeadSele->show();
	connect(_undeadSele, &QPushButton::clicked,
	        this, &Screen::killSelection);

	int bottom = height() - 50;

	_export = new QPushButton("Export all", this);
	_export->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, bottom,
	                     RIGHT_VIEW_WIDTH - 20, 40);
	_export->show();
	connect(_export, &QPushButton::clicked,
	        _list, &ClusterList::exportAll);
	
	bottom -= 50;

	_images = new QPushButton("Export images", this);
	_images->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, bottom,
	                     RIGHT_VIEW_WIDTH - 20, 40);
	_images->show();
	connect(_images, &QPushButton::clicked,
	        this, &Screen::saveImages);

	bottom -= 50;

	_bin.push_back((QWidget **)&_export);
	_bin.push_back((QWidget **)&_invertSele);
	_bin.push_back((QWidget **)&_markSele);
	_bin.push_back((QWidget **)&_unmarkSele);
	_bin.push_back((QWidget **)&_deadSele);
	_bin.push_back((QWidget **)&_undeadSele);
	_bin.push_back((QWidget **)&_newSel);
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

void Screen::markSelection()
{
	bool mark = (QObject::sender() == _markSele);
	Averager *ave = _keeper->getAverager();
	ave->setMarked(mark);
	refreshSelection();
}

void Screen::killSelection()
{
	bool dead = (QObject::sender() == _deadSele);
	Averager *ave = _keeper->getAverager();
	ave->setDead(dead);
	refreshSelection();
}

void Screen::newSelection()
{
	std::vector<MtzFFTPtr> mtzs = _keeper->getMtzsFromSelection();
	_list->makeGroup(mtzs, true);
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
	
	if (_keeper)
	{
		_keeper->getPoints()->repopulate();
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
	if (index == 1 && _selection != NULL)
	{
		_selection->setFocus();
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
