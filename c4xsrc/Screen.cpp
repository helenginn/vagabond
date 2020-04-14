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
#include <QIcon>
#include <iostream>

#include "Screen.h"
#include "SelectionWindow.h"
#include "ClusterList.h"
#include "CorrelLabel.h"
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

	_tabs = NULL;
	_correlLabel = NULL;
	_correlImage = NULL;
	_graph = NULL;
	_keeper = NULL;
	_scroll = NULL;
	_newSel = NULL;
	_markSele = NULL;
	_unmarkSele = NULL;
	_invertSele = NULL;

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
			/*
			_correlLabel->setGeometry(0, tabtop, _tabs->width(), 
			                          _tabs->height());
			*/

			relinkPixmap();
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
				                     _markSele->y(), RIGHT_VIEW_WIDTH - 20, 40);
				_unmarkSele->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
				                     _unmarkSele->y(), RIGHT_VIEW_WIDTH - 20, 40);
				_invertSele->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
				                     _invertSele->y(), RIGHT_VIEW_WIDTH - 20, 40);
				_export->setGeometry(width() - RIGHT_VIEW_WIDTH + 10,
				                      h - 50, RIGHT_VIEW_WIDTH - 20, 40);
			}
		}
	}
}

void Screen::addToolBar()
{
	_toolBar = new QToolBar(this);
	_toolBar->show();
	
	QPushButton *a = new QPushButton("Set average");
	
	QMenu *m = new QMenu(a);
	QAction *a1 = m->addAction("Recalculate");
	connect(a1, &QAction::triggered, this, &Screen::averageGroup);
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
	
	_list->average(static_cast<Averager *>(item));
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
		delete *_bin[i];
		*_bin[i] = 0;
	}

	_bin.clear();
	delete _tabs;
	_tabs = NULL;
}

void Screen::displayResults(Averager *ave)
{
	binTab();
	_tabs = new QTabWidget(this);
	resizeEvent(NULL);

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

	_graph = new QWidget(this);
	_tabs->addTab(_graph, "Axis explorer");

	KeeperGL *gl = new KeeperGL(_graph);
	gl->addAxes();
	gl->setAverager(ave);
	GLPoint *points = gl->getPoints();

	_keeper = gl;

	connect(points, &GLPoint::updateSelection,
	        this, &Screen::refreshSelection);

	_selection = new SelectionWindow(_graph, gl);
	_selection->setPoints(_keeper->getPoints());
	_selection->show();
	_selection->setFocusPolicy(Qt::StrongFocus);

	_scroll = new AxisScroll(_graph);
	_scroll->setAverager(ave);
	_scroll->setPoints(gl->getPoints());
	_scroll->makeLayout();

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


	int bottom = height() - 50;

	_export = new QPushButton("Export all", this);
	_export->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, bottom,
	                     RIGHT_VIEW_WIDTH - 20, 40);
	_export->show();
	connect(_export, &QPushButton::clicked,
	        _list, &ClusterList::exportAll);

	bottom -= 50;

	_bin.push_back((QWidget **)&_export);
	_bin.push_back((QWidget **)&_invertSele);
	_bin.push_back((QWidget **)&_markSele);
	_bin.push_back((QWidget **)&_unmarkSele);
	_bin.push_back((QWidget **)&_scroll);
	_bin.push_back((QWidget **)&_correlLabel);
	_bin.push_back((QWidget **)&_selection);
	_bin.push_back((QWidget **)&_keeper);
	_bin.push_back((QWidget **)&_graph);
	_bin.push_back((QWidget **)&_newSel);

	connect(_tabs, &QTabWidget::tabBarClicked, 
	        this, &Screen::refocus);
	_tabs->show();
	relinkPixmap();
	resizeEvent(NULL);
	refocus(0);
}

void Screen::markSelection()
{
	bool mark = (QObject::sender() == _markSele);
	Averager *ave = _keeper->getAverager();
	ave->setMarked(mark);
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
	
	if (_keeper)
	{
		_keeper->getPoints()->repopulate();
	}
}

void Screen::refocus(int index)
{
	if (index == 1)
	{
		_selection->setFocus();
	}
	else if (index == 0)
	{
		_correlLabel->setFocus();
	}
}

