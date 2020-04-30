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

#include <QMainWindow>
#include <libsrc/mat4x4.h>
#include <libsrc/FFT.h>
#include "MtzFFTPtr.h"

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
class Averager;
class AxisScroll;
class MatrixView;
class KeeperGL;

class Screen : public QMainWindow
{
Q_OBJECT
public:
	Screen(QWidget *widget = NULL);
	
	ClusterList *getList()
	{
		return _list;
	}

	virtual void resizeEvent(QResizeEvent *e);
	void displayResults(Averager *ave);
	void displaySingle(MtzFFTPtr fft);
	
public slots:
	void refocus(int index);
	void changeIndex(int index);
	void saveImages();
	void averageGroup();
	void clusterGroup();
	void newSelection();
	void markSelection();
	void removeCluster();
	void refreshSelection();
	void killSelection();
private:
	void relinkPixmap();
	void addToolBar();
	void binTab();
	void addCorrelImage(Averager *ave);
	void addAxisExplorer(Averager *ave);
	void addCAlphaView();
	void addHKLView(VagFFTPtr fft, std::string filename = "");

	int _currIndex;
	double _scale;
	std::vector<QWidget **> _bin;
	QGraphicsScene *_scene;
	SelectionWindow *_selection;
	MatrixView *_correlImage;
	KeeperGL *_keeper;
	KeeperGL *_hklKeeper;
	KeeperGL *_cAlphaKeeper;
	CorrelLabel *_correlLabel;
	QTabWidget *_tabs;
	QWidget *_graph;
	QWidget *_hkl;
	QWidget *_cAlpha;
	AxisScroll *_scroll;
	QTreeWidget *_inputTree;
	QToolBar *_toolBar;
	QPushButton *_newSel;
	QPushButton *_markSele;
	QPushButton *_unmarkSele;
	QPushButton *_deadSele;
	QPushButton *_undeadSele;
	QPushButton *_invertSele;
	QPushButton *_export;
	QPushButton *_images;
	QPushButton *_toggleDead;
	QPlainTextEdit *_ucLabel;
	QAction *_cluster;
	ClusterList *_list;

	mat4x4 _storeHKL;
	mat4x4 _storeCAlpha;
};
