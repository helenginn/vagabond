//
//  VagWindow.cpp
//  Vagabond
//
//  Created by Helen Ginn on 11/12/2016.
//  Copyright (c) 2016 Helen Ginn. All rights reserved.
//

#define DEFAULT_WIDTH 900
#define DEFAULT_HEIGHT 700

#include "../../libsrc/shared_ptrs.h"
#include "VagWindow.h"
#include <QtCore/qdebug.h>
#include <QtCore/qalgorithms.h>
#include <QtWidgets/qgraphicsitem.h>
#include <QtWidgets/qmenubar.h>
#include <QtWidgets/qmessagebox.h>
#include <iostream>
#include "../../libsrc/Options.h"
#include "InstructionThread.h"
#include "Dialogue.h"
#include "../../libsrc/Crystal.h"

void VagWindow::makeButtons()
{
	bSuperimpose = new QPushButton("Superimpose", this);
	bSuperimpose->setGeometry(700, 0, 200, 50);
	bSuperimpose->setEnabled(false);
	connect(bSuperimpose, SIGNAL(clicked()), this, SLOT(pushSuperimpose()));
	
	bRefinePos = new QPushButton("Refine positions to PDB", this);
	bRefinePos->setGeometry(700, 50, 200, 50);
	bRefinePos->setEnabled(false);
	connect(bRefinePos, SIGNAL(clicked()), this, SLOT(pushRefinePositions()));

	bRefineFlex = new QPushButton("Refine flexibility to PDB", this);
	bRefineFlex->setGeometry(700, 100, 200, 50);
	bRefineFlex->setEnabled(false);
	connect(bRefineFlex, SIGNAL(clicked()), this, SLOT(pushRefineFlexibility()));

	bChangeBMult = new QPushButton("Set hetatm B multiplier", this);
	bChangeBMult->setGeometry(700, 200, 200, 50);
	bChangeBMult->setEnabled(false);
	connect(bChangeBMult, SIGNAL(clicked()), this, SLOT(pushBMultiplier()));
	
	bExploreMolecule = new QPushButton("Explore (first) molecule", this);
	bExploreMolecule->setGeometry(700, 300, 200, 50);
	bExploreMolecule->setEnabled(false);
	connect(bExploreMolecule, SIGNAL(clicked()), this, SLOT(pushExploreMcule()));
	
	_myDialogue = NULL;
	_explorer = NULL;
}

VagWindow::VagWindow(QWidget *parent,
                     int argc, char *argv[]) : QMainWindow(parent)
{
	this->resize(DEFAULT_WIDTH, DEFAULT_HEIGHT);
    
    this->setWindowTitle("Vagabond");
	
	display = new VagabondGLWidget(this);
	display->setGeometry(0, 0, 700, 700);
	
	makeButtons();
	
	_argc = argc;
	_argv = argv;
	
	_instructionThread.setVagWindow(this);
	_instructionThread.start();
	
	_instructionType = InstructionTypeNone;
}

void VagWindow::waitForInstructions()
{
	OptionsPtr options = OptionsPtr(new Options(_argc, (const char **)_argv));
	options->setNotify((Notifiable *)this);
	Options::setRuntimeOptions(options);
	options->run();
	
    while (true)
    {
        mutex.lock();
        wait.wait(&mutex);
        
	    switch (_instructionType)
	    {
	        case InstructionTypeSuperimpose:
	        options->superimposeAll();
	        break;

	        case InstructionTypeRefinePositions:
	        options->refineAll(RefinementModelPos, 1);
	        break;
	        
	        case InstructionTypeRefineFlexibility:
	        options->refineAll(RefinementFlexibility, 1);
	        break;
	    
	        case InstructionTypeChangeBMult:
	        options->applyBMultiplier();
	        break;
	    
	        default:
	        break;
	    }
	    
        mutex.unlock();
    }
}

void VagWindow::enable()
{
    bSuperimpose->setEnabled(true);
    bRefinePos->setEnabled(true);
    bRefineFlex->setEnabled(true);
    bChangeBMult->setEnabled(true);
    bExploreMolecule->setEnabled(true);
}

void VagWindow::disable()
{
    bSuperimpose->setEnabled(false);
    bRefinePos->setEnabled(false);
    bRefineFlex->setEnabled(false);
    bChangeBMult->setEnabled(false);
    bExploreMolecule->setEnabled(false);
}

void VagWindow::resizeEvent(QResizeEvent *event)
{
    int w = width();
    int h = height();
    int max_w = w;
    int max_h = h;
    
    h = (h < w) ? h : w;
    w = (w < h) ? w : h;
    
    int top = max_h / 2 - h / 2;
    int left = max_w / 2 - w / 2;
}

void VagWindow::pushSuperimpose()
{
    _instructionType = InstructionTypeSuperimpose;
	wait.wakeAll();
}

void VagWindow::pushRefinePositions()
{
    _instructionType = InstructionTypeRefinePositions;
	wait.wakeAll();
}

void VagWindow::pushBMultiplier()
{
	delete _myDialogue;
    _myDialogue = new Dialogue(this, "New B multiplier",
    								"Enter scale factor for hetatm Bs",
    								"1.0",
    								"Set B multiplier");
	_myDialogue->setTag(DialogueBMultiplier);
    _myDialogue->setWindow(this);
	_myDialogue->show();
}

void VagWindow::pushRefineFlexibility()
{
    _instructionType = InstructionTypeRefineFlexibility;
	wait.wakeAll();
}

void VagWindow::pushExploreMcule()
{
    delete _explorer;
    _explorer = NULL;
    
    OptionsPtr options = Options::getRuntimeOptions();
    int crystalCount = options->crystalCount();
    
    if (crystalCount)
    {
        int moleculeCount = options->getCrystal(0)->moleculeCount();
        
        if (moleculeCount)
        {
            MoleculePtr molecule = options->getCrystal(0)->molecule(0);
            _explorer = new MoleculeExplorer(this, molecule);
            _explorer->show();
        }
    }
}


void VagWindow::receiveDialogue(DialogueType type, std::string diagString)
{
	std::vector<double> trial;
	
	while (true)
	{
		size_t pos = diagString.find(" ");
		
		if (pos >= diagString.length()) pos = diagString.length();
		
		std::string substr = diagString.substr(0, pos);
		double value = atof(substr.c_str());
		trial.push_back(value);
		
		if (pos >= diagString.length()) break;
		diagString = diagString.substr(pos + 1, diagString.length() - pos - 1);
	}

    if (!trial.size())
    {
        return;
    }
	
	if (type == DialogueBMultiplier)
	{
		OptionsPtr options = Options::getRuntimeOptions();
		options->setBMult(trial[0]);
		_instructionType = InstructionTypeChangeBMult;
        wait.wakeAll();		
	}
	
	_myDialogue->hide();
	delete _myDialogue;
	_myDialogue = NULL;
}


VagWindow::~VagWindow()
{
    delete bSuperimpose;
    delete bRefinePos;
    delete bRefineFlex;
    delete _myDialogue;
    delete _explorer;
}
