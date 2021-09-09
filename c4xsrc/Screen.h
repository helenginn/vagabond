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

#ifndef __c4xsrc__screen__
#define __c4xsrc__screen__

#include <QMainWindow>
#include <hcsrc/mat4x4.h>
#include <FFT.h>
#include "MtzFFTPtr.h"
#include "PlotView.h"

class CorrelLabel;
class MtzFile;
class QPlainTextEdit;
class QTreeWidget;
class QImage;
class QPushButton;
class QGraphicsScene;
class SelectionWindow;
class QLabel;
class ClusterList;
class Group;
class MatrixView;
class KeeperGL;
class PlotView;
class QKeyEvent;
class ColumnView;

typedef struct
{
	QAction *amp;
	QAction *ca;
	QAction *uc;
	QAction *vec;
	QAction *csv;
	QAction *top;
	QAction *orig;
	QAction *self;
} GroupUses;

class C4XAcceptor
{
public:
	virtual void finished() = 0;
	virtual ~C4XAcceptor() {}

};

class Screen : public QMainWindow
{
Q_OBJECT
public:
	Screen(QWidget *widget = NULL);
	
	ClusterList *getList()
	{
		return _list;
	}

	void displayResults(Group *ave);
	void displaySingle(MtzFFTPtr fft);
	void updateToolbar(Group *grp);
	void addCSVSwitcher();
	void displaySettings();
	
	void setReturnJourney(C4XAcceptor *ptr);
	
signals:
	void refreshed();
public slots:
	void refocus(int index);
	void changeIndex(int index);
	void saveImages();
	void writeAverageMTZ();
	void clusterGroup();
	void newSelection();
	void markSelection();
	void removeCluster();
	void rotateDegrees();
	void refreshSelection();
	void collapsePositions();
	void reorder();
	void reindex();
	void coverage();
	void changeColour();
	void killSelection();
	void reorderByFile();
	void returnToSender();
	void exportText();
	void plotSpin();
protected:

	virtual void keyPressEvent(QKeyEvent *e);
	virtual void keyReleaseEvent(QKeyEvent *e);
private:
	void addSideButtons();
	void relinkPixmap();
	void addColour(QString colour, QString display, QMenu *m);
	void addToolBar();
	void binTab();
	void addCorrelImage(Group *ave);
	void addPlotView(PlotView **view, Group *ave,
                         std::string title, PlotType type);

	void addCAlphaView();
	void addColumnView(Group *ave);
	void addHKLView(VagFFTPtr fft, std::string filename = "");

	void addSideButton(QWidget **buttPtr, std::string title, int *top);

	int _currIndex;
	double _scale;
	QMenu *_menu;
	std::vector<QWidget **> _bin;
	QGraphicsScene *_scene;
	MatrixView *_correlImage;
	PlotView *_svdView;
	PlotView *_ucView;
	PlotView *_rView;
	Group *_group;
	KeeperGL *_hklKeeper;
	KeeperGL *_cAlphaKeeper;
	CorrelLabel *_correlLabel;
	QTabWidget *_tabs;
	QWidget *_hkl;
	QWidget *_cAlpha;
	QWidget *_side;
	QTreeWidget *_inputTree;
	QToolBar *_toolBar;
	QPushButton *_newSel;
	QPushButton *_invertSele;
	QPushButton *_markSele;
	QPushButton *_deadSele;
	QPushButton *_changeColour;
	QPushButton *_export;
	QPushButton *_images;
	QPushButton *_toggleDead;
	QPushButton *_changeData;
	QPushButton *_coverage;
	QPushButton *_reorder;
	QPlainTextEdit *_ucLabel;
	QAction *_cluster;
	ClusterList *_list;
	ColumnView *_columnView;
	
	GroupUses uses;
	
	C4XAcceptor *_returnJourney;

	mat4x4 _storeHKL;
};


#endif
