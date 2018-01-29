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
    InstructionTypeOpenPDB,
    InstructionTypeOpenMTZ,
    InstructionTypeSuperimpose,
    InstructionTypeRefinePositions,
    InstructionTypeRefineFlexibility,
    InstructionTypeChangeBMult,
    InstructionTypeRecalculateFFT,
    InstructionTypeSetOutputDir,
} InstructionType;


class VagWindow : public QMainWindow, public Notifiable
{
    Q_OBJECT
    
public:
    VagWindow(QWidget *parent = 0, int argc = 0, char *argv[] = NULL);
    ~VagWindow();

    virtual void disable();
    virtual void enable();
    void waitForInstructions();
    void receiveDialogue(DialogueType type, std::string diagString);

    virtual void setMessage(std::string message);

protected:
    virtual void resizeEvent(QResizeEvent *event);
private slots:
    void pushSuperimpose();
    void pushRefinePositions();
    void pushRefineFlexibility();
    void pushBMultiplier();
    void pushExploreMcule();
    void recalculateFFT();

    void openPDB();
    void openMTZ();
    void setOutput();

private:
    VagabondGLWidget *display;
    QWaitCondition wait;
    QMutex mutex;
    InstructionThread _instructionThread;
    InstructionType _instructionType;
    Dialogue *_myDialogue;
    MoleculeExplorer *_explorer;
    QFileDialog *_fileDialogue;   
    
    void updateExplorerButton();
    QLabel *_lStatus;

    /* Buttons down the side */
    QPushButton *bSuperimpose;
    QPushButton *bRefinePos;
    QPushButton *bRefineFlex;
    QPushButton *bChangeBMult;
    QPushButton *bExploreMolecule;
    QPushButton *bRecalculate;

    std::string _outputDir;
    std::string _pdbName; 
    std::string _mtzName;
    int _argc;
    char **_argv;

    void makeMenu();
    void makeButtons();
};

#endif /* defined(__Vagabond__VagWindow__) */
