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
#include "VagabondGLWidget.h"
#include "../../libsrc/Notifiable.h"
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

typedef enum
{
    InstructionTypeNone,
    InstructionTypeSuperimpose,
    InstructionTypeRefinePositions,
    InstructionTypeRefineFlexibility,
    InstructionTypeChangeBMult,
} InstructionType;


class VagWindow : public QMainWindow, public Notifiable
{
    Q_OBJECT
    
public:
    VagWindow(QWidget *parent = 0, int argc = 0, char *argv[] = NULL);

	/* Buttons down the side */
    QPushButton *bSuperimpose;
    QPushButton *bRefinePos;
    QPushButton *bRefineFlex;
    QPushButton *bChangeBMult;
    QPushButton *bExploreMolecule;

    ~VagWindow();
    
    void makeButtons();
    virtual void disable();
    virtual void enable();
    void waitForInstructions();
    void receiveDialogue(DialogueType type, std::string diagString);
    
protected:
    virtual void resizeEvent(QResizeEvent *event);
private slots:
	void pushSuperimpose();
	void pushRefinePositions();
	void pushRefineFlexibility();
	void pushBMultiplier();
	void pushExploreMcule();
	
private:
    VagabondGLWidget *display;
    QWaitCondition wait;
    QMutex mutex;
    InstructionThread _instructionThread;
    InstructionType _instructionType;
    Dialogue *_myDialogue;
    MoleculeExplorer *_explorer;
    
    int _argc;
    char **_argv;
};

#endif /* defined(__Vagabond__VagWindow__) */
