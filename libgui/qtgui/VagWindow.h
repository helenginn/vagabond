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

#include "../../libsrc/Shouter.h"
#include "InstructionThread.h"
#include "Dialogue.h"
#include "StartScreen.h"
#include "CrystalExplorer.h"
#include "ErroneousZone.h"
#include "../../libsrc/Notifiable.h"

class VagabondGLWidget;
#include "VagabondGLWidget.h"

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

	void gotoResidueDialogue();
	void receiveGotoResidue(std::string diagString);

	virtual void setMessage(std::string message);
	virtual void wakeup();
	virtual void setRenderDensity();
	virtual void focusOnPosition(vec3 pos);
	virtual void pause(bool on);
	virtual void appendToLog(char *msg);
	void pushSendInstruction(InstructionType inst);
	
	void toggleLog();
	
	void setStartScreen(StartScreen *screen)
	{
		_startScreen = screen;
	}
	
	InstructionThread *getInstructionThread()
	{
		return &_instructionThread;	
	}
	
	void setInstructionType(InstructionType type)
	{
		_instructionType = type;
	}
	
	VagabondGLWidget *getDisplay()
	{
		return display;
	}
	
	void setObject(void *obj)
	{
		_obj = obj;
	}

	int attemptLoadAndRun();
	
signals:
	void errorReceived(Shouter *shout);
	void appendSignal();
protected:
	virtual void resizeEvent(QResizeEvent *);
private slots:
	void pushBMultiplier();
	void pushExploreCrystal();
	void restorePreviousState();

	void toggleScaling(ScalingType type);
	void adjustBFactor();
	void refitBackbone();
	void fixErroneousZones();
	void append();
	void openInCoot();
	void addPDBFile();

	void displayMessage(Shouter *shout);
private:
	VagabondGLWidget *display;
	void fixLabelChoice(LabelChoice &choice);
	QWaitCondition wait;
	QMutex mutex;
	InstructionThread _instructionThread;
	Dialogue *_myDialogue;
	CrystalExplorer *_xtalExplorer;
	ErroneousZone *_errorExplorer;
	QFileDialog *_fileDialogue;   
	QTextEdit *_logView;

	void dialogueModify(Setter set, std::string title, 
	                    double _default = 100);
	void splitBond();
	void resetSides();
	void resetMotion();
	void refineToEnd();
	void modelPosToEnd();
	void sidechainsToEnd();
	void fitWholeMolecules();
	void getPolymerMonomerCrystal(PolymerPtr *poly, CrystalPtr *cryst, 
	                              MonomerPtr *monomer);
	void menuItem(QMenu *menu, std::string title,
                         InstructionType instr);

	QLabel *_lStatus;
	StartScreen *_startScreen;

	/* Buttons down the side */
	QPushButton *bRefinePos;
	QPushButton *bChelate;
	QPushButton *bChangeBMult;
	QPushButton *bExploreCrystal;
	QPushButton *bRecalculate;
	QPushButton *bRefineDensity;
	QPushButton *bFitWholeR;
	QPushButton *bFitWholeT;
	QPushButton *bPrevious;
	QPushButton *bCoot;
	std::vector<QPushButton *> buttons;
	
	/* Scaling menu options */
	QAction *_qaShell, *_qaKB, *_qaK;
	
	InstructionThread *_thread;
	
	std::vector<QMenu *> menus;
	std::vector<QAction *> actions;

	void displayScaling();
	
	int _rangeStart;
	int _rangeEnd;

	bool _showingLog;
	int _argc;
	char **_argv;
	
	std::mutex _guiOutMut;
	std::string _guiOut;
	std::string _pdbStr;

	void makeMenu();
	void makeButtons();
	void *_obj;
};

#endif /* defined(__Vagabond__VagWindow__) */
