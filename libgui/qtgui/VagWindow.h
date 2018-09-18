//
//  VagWindow.h
//  Vagabond
//  I am completely aware of the filename hilarity
//
//  Created by Helen Ginn on 11/12/2016.
//  Copyright (c) 2017 Helen Ginn. All rights reserved.
//

#ifndef __Vagabond__VagWindow__
#define __Vagabond__VagWindow__

#include <stdio.h>
#include <QtCore/qglobal.h>
#include <QtWidgets/qapplication.h>
#include <QtWidgets/qmainwindow.h>
#include <QtWidgets/qwidget.h>
#include <QtWidgets/qpushbutton.h>
#include <QtWidgets/qlabel.h>
#include <QtGui/qpixmap.h>
#include <QtWidgets/qfiledialog.h>
#include <QtWidgets/qgraphicsview.h>
#include <vector>
#include <QtCore/qsignalmapper.h>
#include <QtCore/qwaitcondition.h>
#include <QtCore/qmutex.h>

#include "InstructionThread.h"
#include "Dialogue.h"
#include "MoleculeExplorer.h"
#include "CrystalExplorer.h"
#include "VagabondGLWidget.h"
#include "../../libsrc/Notifiable.h"


class VagWindow : public QMainWindow, public Notifiable
{
	Q_OBJECT

public:
	VagWindow(QWidget *parent = 0, int argc = 0, char *argv[] = NULL);
	~VagWindow();

	virtual void disable();
	virtual void enable();
	int waitForInstructions();
	virtual bool isRunningSomething();
	void receiveDialogue(DialogueType type, std::string diagString);

	virtual void setMessage(std::string message);
	virtual void wakeup();
	virtual void setRenderDensity();
	
	InstructionThread *getInstructionThread()
	{
		return &_instructionThread;	
	}
	
protected:
	virtual void resizeEvent(QResizeEvent *);
private slots:
	void pushSuperimpose();
	void pushFitWholeR();
	void pushFitWholeT();
	void pushRefinePositions();
	void pushRefineFlexibility();
	void pushBMultiplier();
	void pushExploreMcule();
	void pushExploreCrystal();
	void pushRefineDensity();
	void recalculateFFT();
	void openInCoot();
	void findDisulphides();
	void pushBackboneAnalysis();
	void restorePreviousState();
	void refineWaterNetwork();

	void toggleScaling(ScalingType type);
	void adjustBFactor();

private:
	VagabondGLWidget *display;
	QWaitCondition wait;
	QMutex mutex;
	InstructionThread _instructionThread;
	Dialogue *_myDialogue;
	MoleculeExplorer *_moleExplorer;
	CrystalExplorer *_xtalExplorer;
	QFileDialog *_fileDialogue;   

	void splitBond();
	void squeezeToEnd();
	void updateExplorerButton();
	void refineToEnd();
	void modelPosToEnd();
	void sidechainsToEnd();
	void getPolymerMonomerCrystal(PolymerPtr *poly, CrystalPtr *cryst, MonomerPtr *monomer);
	QLabel *_lStatus;

	/* Buttons down the side */
	QPushButton *bSuperimpose;
	QPushButton *bRefinePos;
	QPushButton *bBackbone;
	QPushButton *bChangeBMult;
	QPushButton *bExploreMolecule;
	QPushButton *bExploreCrystal;
	QPushButton *bRecalculate;
	QPushButton *bRefineDensity;
	QPushButton *bFitWholeR;
	QPushButton *bFitWholeT;
	QPushButton *bWaterNetwork;
	QPushButton *bPrevious;
	QPushButton *bCoot;
	std::vector<QPushButton *> buttons;
	
	/* Scaling menu options */
	QAction *_qaShell, *_qaKB, *_qaK;
	std::vector<QMenu *> menus;
	std::vector<QAction *> actions;

	void displayScaling();

	int _argc;
	char **_argv;

	void makeMenu();
	void makeButtons();
};

#endif /* defined(__Vagabond__VagWindow__) */
