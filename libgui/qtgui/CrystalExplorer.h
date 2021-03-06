//  Vagabond
//  CrystalExplorer.h
//
//  Created by Helen Ginn on 4/8/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#ifndef __Vagabond__CrystalExplorer__
#define __Vagabond__CrystalExplorer__

#include <stdio.h>

#include "../../libsrc/Notifiable.h"
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
class VagabondGLWidget;
class VagWindow;

class CrystalExplorer : public QMainWindow
{
	Q_OBJECT
	
public:
	CrystalExplorer(QWidget *parent = 0,
	                CrystalPtr crystal = CrystalPtr());
	
	~CrystalExplorer();
	
	void setKeeper(VagabondGLWidget *keeper)
	{
		_keeper = keeper;
	}
	
	void setVagWindow(VagWindow *window)
	{
		_vagWindow = window;
	}
	
	void updateCorrelation();
private slots:
	void clickedMoleListItem();
	void slideScale(QSlider *s);
	void pushWritePNG();
	void pushSequence();
	void pushFitSides();
	void pushFitIntra();
	void pushBackbone();
	void refineSponges();
	void pushFitMotion();
	void pushRigidBody();
	void pushResetSides();
	void pushResetIntra();
	void pushResetMotion();
	void pushDistanceMatrix();

private:
	void sendObject(InstructionType type);

	MoleculePtr _currMole;
	MoleculeExplorer *_moleExplorer;
	
	VagabondGLWidget *_keeper;

	CrystalPtr _crystal;
	QListWidget *_moleList;
	VagWindow *_vagWindow;
	
	std::vector<QWidget *> _widgets;
	void populateList();
	void clearWidgets();
};

#endif
