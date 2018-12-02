//  Vagabond
//  CrystalExplorer.h
//
//  Created by Helen Ginn on 4/8/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#ifndef __Vagabond__CrystalExplorer__
#define __Vagabond__CrystalExplorer__

#include <stdio.h>

#include <QtCore/qglobal.h>
#include <QtWidgets/qapplication.h>
#include <QtWidgets/qwidget.h>
#include <QtWidgets/qpushbutton.h>
#include <QtWidgets/qlabel.h>
#include <QtWidgets/qtextedit.h>
#include <QtWidgets/qmainwindow.h>
#include <QtWidgets/qlistwidget.h>
#include "../../libsrc/shared_ptrs.h"

class MoleculeExplorer;
class GLKeeper;

class CrystalExplorer : public QMainWindow
{
	Q_OBJECT
	
public:
	CrystalExplorer(QWidget *parent = 0,
	                CrystalPtr crystal = CrystalPtr());
	
	~CrystalExplorer();
	
	void setKeeper(GLKeeper *keeper)
	{
		_keeper = keeper;
	}
private slots:
	void clickedMoleListItem();
	void pushSequence();
	
private:
	MoleculePtr _currMole;
	MoleculeExplorer *_moleExplorer;
	
	GLKeeper *_keeper;

	CrystalPtr _crystal;
	QListWidget *_moleList;
	std::vector<QWidget *> _widgets;
	void populateList();
	void clearWidgets();
};

#endif
