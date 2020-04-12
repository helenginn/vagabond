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
#include <QPushButton>
#include <QMessageBox>
#include <QMenu>
#include <QLabel>
#include <QIcon>
#include <iostream>

#include "Screen.h"
#include "SelectionWindow.h"
#include "ClusterList.h"
#include "KeeperGL.h"
#include "GLPoint.h"
#include "Averager.h"
#include "AxisScroll.h"
#include "MatrixView.h"
#include "MtzFFTPtr.h"

#define TOOL_BAR_HEIGHT 50
#define TREE_VIEW_WIDTH 200
#define TAB_VIEW_WIDTH 500
#define RIGHT_VIEW_WIDTH 200

Screen::Screen(QWidget *widget) : QMainWindow(widget)
{
	setGeometry(0, 0, 1000, 600);

	_tabs = NULL;
	_correlLabel = NULL;
	_graph = NULL;
	_keeper = NULL;
	_scroll = NULL;

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
			_correlLabel->setGeometry(0, 0, _tabs->width(), _tabs->height());
			QImage i = _correlImage->scaled(_correlLabel->width(), 
			                                _correlLabel->height(),
			                                Qt::KeepAspectRatioByExpanding);
			_correlLabel->setPixmap(QPixmap::fromImage(i));
		}
		
		if (_graph)
		{
			_graph->setGeometry(0, 0, _tabs->width(), _tabs->height());
			
			if (_keeper)
			{
				_keeper->setGeometry(0, 60, _tabs->width(), 
				                     _tabs->height() - 60);
			}
			
			if (_selection)
			{
				_selection->setGeometry(-2, 58, _tabs->width() + 2, 
				                        _tabs->height() - 60);
			}
			
			if (_scroll)
			{
				_scroll->setGeometry(0, 0, _tabs->width(), 60);
			}

			if (_newSel)
			{
				_newSel->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, 10,
				                     RIGHT_VIEW_WIDTH - 20, 40);
			}
		}
	}
}

void Screen::addToolBar()
{
	_toolBar = new QToolBar(this);
	_toolBar->show();
	
	QAction *a = _toolBar->addAction("Average");
	connect(a, &QAction::triggered, this, &Screen::averageGroup);

	QAction *c = _toolBar->addAction("Cluster");
	connect(c, &QAction::triggered, this, &Screen::clusterGroup);

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

	if (ave->getCorrelMatrix())
	{
		_correlImage = ave->getCorrelMatrix();
		_correlImage->updateSelection();

		QLabel *l = new QLabel(NULL);
		l->setGeometry(0, 0, _tabs->height(), _tabs->width());
		QImage i = _correlImage->scaled(l->width() - 40, l->height(),
		                                Qt::KeepAspectRatioByExpanding);
		l->setPixmap(QPixmap::fromImage(i));
		_tabs->addTab(l, "Correlation matrix");

		_correlLabel = l;

		_graph = new QWidget(this);
		_graph->setGeometry(0, 0, _tabs->width(), _tabs->width());
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
		
		_newSel = new QPushButton("New selection", this);
		_newSel->setGeometry(width() - RIGHT_VIEW_WIDTH + 10, top,
		                     RIGHT_VIEW_WIDTH - 20, 40);
		_newSel->show();
		
		top += 50;
		
		QMenu *menu = new QMenu(this);
		_withAve = menu->addAction("Using existing average");
		connect(_withAve, &QAction::triggered,
		        this, &Screen::newSelection);
		_newAve = menu->addAction("Calculate new average");
		connect(_newAve, &QAction::triggered,
		        this, &Screen::newSelection);
		_newSel->setMenu(menu);

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

		_bin.push_back((QWidget **)&_invertSele);
		_bin.push_back((QWidget **)&_markSele);
		_bin.push_back((QWidget **)&_unmarkSele);
		_bin.push_back((QWidget **)&_newAve);
		_bin.push_back((QWidget **)&_withAve);
		_bin.push_back((QWidget **)&_scroll);
		_bin.push_back((QWidget **)&_correlLabel);
		_bin.push_back((QWidget **)&_selection);
		_bin.push_back((QWidget **)&_keeper);
		_bin.push_back((QWidget **)&_graph);
		_bin.push_back((QWidget **)&_newSel);
	}

	connect(_tabs, &QTabWidget::tabBarClicked, 
	        this, &Screen::refocus);
	_tabs->show();
}

void Screen::markSelection()
{
	bool mark = (QObject::sender() == _markSele);
	std::vector<MtzFFTPtr> mtzs = _keeper->getMtzs();
	_list->markMtzs(mtzs, mark);
	
	Averager *ave = _keeper->getAverager();
	QFont curr = ave->font(0);
	curr.setBold(mark);
	ave->setFont(0, curr);
}

void Screen::newSelection()
{
	bool withAve = (QObject::sender() == _withAve);
	std::vector<MtzFFTPtr> mtzs = _keeper->getMtzsFromSelection();
	_list->makeGroup(mtzs, withAve);
}

void Screen::refreshSelection()
{
	if (_correlImage)
	{
		_correlImage->updateSelection();
		_correlLabel->setPixmap(QPixmap::fromImage(*_correlImage));
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
		_keeper->setFocus();
	}
}
