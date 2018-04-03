//
//  VagWindow.cpp
//  Vagabond
//
//  Created by Helen Ginn on 11/12/2016.
//  Copyright (c) 2016 Helen Ginn. All rights reserved.
//

#define DEFAULT_WIDTH 900
#define DEFAULT_HEIGHT 660
#define STATUS_HEIGHT 60
#define BUTTON_WIDTH 200

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
#include "../../libsrc/Monomer.h"
#include "../../libsrc/Polymer.h"
#include "ChainMenuAction.h"
#include "../../libsrc/FileReader.h"

void VagWindow::makeMenu()
{
    _fileDialogue = NULL;

    QMenu *fileMenu = menuBar()->addMenu(tr("&File"));
    QAction *aOutput = fileMenu->addAction(tr("Set &output"));
    connect(aOutput, &QAction::triggered, this, &VagWindow::setOutput);
    QAction *aOpenPDB = fileMenu->addAction(tr("Open &model..."));
    connect(aOpenPDB, &QAction::triggered, this, &VagWindow::openModel);
    QAction *aOpenMTZ = fileMenu->addAction(tr("Open &MTZ..."));
    connect(aOpenMTZ, &QAction::triggered, this, &VagWindow::openMTZ);
}

void VagWindow::makeButtons()
{
	buttons.clear();

    bSuperimpose = new QPushButton("Superimpose", this);
    bSuperimpose->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 0, BUTTON_WIDTH , 50);
    bSuperimpose->setEnabled(false);
    connect(bSuperimpose, SIGNAL(clicked()), this, SLOT(pushSuperimpose()));
	buttons.push_back(bSuperimpose);
    
    bFitWholeT = new QPushButton("Fit molecule translation", this);
    bFitWholeT->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 50, BUTTON_WIDTH , 25);
    bFitWholeT->setEnabled(false);
    connect(bFitWholeT, SIGNAL(clicked()), this, SLOT(pushFitWholeT()));
	buttons.push_back(bFitWholeT);
    
    bFitWholeR = new QPushButton("Fit molecule rotation", this);
    bFitWholeR->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 75, BUTTON_WIDTH , 25);
    bFitWholeR->setEnabled(false);
    connect(bFitWholeR, SIGNAL(clicked()), this, SLOT(pushFitWholeR()));
	buttons.push_back(bFitWholeR);
    
    bRefinePos = new QPushButton("Refine positions to PDB", this);
    bRefinePos->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 100, BUTTON_WIDTH , 50);
    bRefinePos->setEnabled(false);
    connect(bRefinePos, SIGNAL(clicked()), this, SLOT(pushRefinePositions()));
	buttons.push_back(bRefinePos);

    bBackbone = new QPushButton("Refine backbone", this);
    bBackbone->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 150, BUTTON_WIDTH , 50);
    bBackbone->setEnabled(false);
    connect(bBackbone, SIGNAL(clicked()), this, SLOT(pushBackboneAnalysis()));
	buttons.push_back(bBackbone);

    bRefineDensity = new QPushButton("Refine sidechains to density", this);
    bRefineDensity->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 200, BUTTON_WIDTH , 50);
    bRefineDensity->setEnabled(false);
    connect(bRefineDensity, SIGNAL(clicked()), this, SLOT(pushRefineDensity()));
	buttons.push_back(bRefineDensity);

    bRecalculate = new QPushButton("Recalculate FFT", this);
    bRecalculate->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 250, BUTTON_WIDTH , 50);
    bRecalculate->setEnabled(false);
    connect(bRecalculate, SIGNAL(clicked()), this, SLOT(recalculateFFT()));    
	buttons.push_back(bRecalculate);

    bChangeBMult = new QPushButton("Set hetatm B multiplier", this);
    bChangeBMult->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 300, BUTTON_WIDTH , 50);
    bChangeBMult->setEnabled(false);
    connect(bChangeBMult, SIGNAL(clicked()), this, SLOT(pushBMultiplier()));
	buttons.push_back(bChangeBMult);
    
    bFindSS = new QPushButton("Find disulphides", this);
    bFindSS->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 350, BUTTON_WIDTH , 50);
    bFindSS->setEnabled(false);
    connect(bFindSS, SIGNAL(clicked()), this, SLOT(findDisulphides()));
	buttons.push_back(bFindSS);
    
    bExploreMolecule = new QPushButton("Explore molecule", this);
    bExploreMolecule->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 400, BUTTON_WIDTH , 50);
    bExploreMolecule->setEnabled(false);
    bExploreMolecule->setMenu(new QMenu(this));
	buttons.push_back(bExploreMolecule);
    
    bPrevious = new QPushButton("Restore previous state", this);
    bPrevious->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 450, BUTTON_WIDTH , 50);
    bPrevious->setEnabled(false);
    connect(bPrevious, SIGNAL(clicked()), this, SLOT(restorePreviousState()));
	buttons.push_back(bPrevious);
    
    _myDialogue = NULL;
    _explorer = NULL;
}

void VagWindow::updateExplorerButton()
{
    OptionsPtr options = Options::getRuntimeOptions();
    
    if (!options->crystalCount())
    {
        return;
    }

    CrystalPtr crystal = options->getActiveCrystal();

    QMenu *moleculeMenu = new QMenu(bExploreMolecule);
    
    for (int i = 0; i < crystal->moleculeCount(); i++)
    {
        MoleculePtr molecule = crystal->molecule(i);
        
        if (!molecule->isPolymer())
        {
            continue;
        }

        ChainMenuAction *action = new ChainMenuAction(moleculeMenu, molecule);
        connect(action, &QAction::triggered, this, &VagWindow::pushExploreMcule);
    }

    bExploreMolecule->setMenu(moleculeMenu);
    connect(bExploreMolecule, SIGNAL(clicked()), bExploreMolecule, SLOT(showMenu())); 
}

void VagWindow::resizeEvent(QResizeEvent *)
{
    int w = width();
    int h = height();

	display->setGeometry(0, 0, w - BUTTON_WIDTH, h - STATUS_HEIGHT);
	display->resizeGL();
	_lStatus->setGeometry(10, h - STATUS_HEIGHT, w - BUTTON_WIDTH, STATUS_HEIGHT);

	for (int i = 0; i < buttons.size(); i++)
	{
		buttons[i]->move(w - BUTTON_WIDTH, buttons[i]->y());
	}
}

VagWindow::VagWindow(QWidget *parent,
                     int argc, char *argv[]) : QMainWindow(parent)
{
    this->resize(DEFAULT_WIDTH, DEFAULT_HEIGHT);
    
    this->setWindowTitle("Vagabond");
    
    display = new VagabondGLWidget(this);
	display->setGeometry(0, 0, DEFAULT_WIDTH - BUTTON_WIDTH, DEFAULT_HEIGHT -
	                     STATUS_HEIGHT);
    
    makeButtons();
    makeMenu();    

    QFont bigFont = QFont("Helvetica", 16);
    _lStatus = new QLabel("Vagabond at your service.", this);
    _lStatus->setFont(bigFont);
    _lStatus->setGeometry(10, DEFAULT_HEIGHT - STATUS_HEIGHT, 680, STATUS_HEIGHT);
    _lStatus->show();

    _argc = argc;
    _argv = argv;
    
    _instructionType = InstructionTypeNone;
    _instructionThread.setVagWindow(this);
    _instructionThread.start(); 

	display->setFocus();
	display->setFocusPolicy(Qt::StrongFocus);
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

			case InstructionTypeBackboneAnalysis:
			options->backboneAnalysis();
			break;

            case InstructionTypeFitWholeMoleculeTranslation: 
            options->fitWholeMolecule(true, false);
            break;

            case InstructionTypeFitWholeMoleculeRotation: 
            options->fitWholeMolecule(false, true);
            break;

            case InstructionTypeSqueezeToEnd: 
            squeezeToEnd();
            break;

            case InstructionTypeSidechainsToEnd: 
            sidechainsToEnd();
            break;

            case InstructionTypeModelPosToEnd: 
            modelPosToEnd();
            break;

            case InstructionTypeRefineToEnd: 
            refineToEnd();
            break;

            case InstructionTypeRefineDensity: 
            options->refineAll(RefinementSidechain, 1);
            break;

            case InstructionTypeChangeBMult:
            options->applyBMultiplier();
            break;
       
            case InstructionTypeOpenPDB:
            options->openModel(_pdbName);
            updateExplorerButton();
            break;

            case InstructionTypeOpenMTZ:
            options->openMTZ(_mtzName);
            break;

            case InstructionTypeFindDisulphides:
            options->findDisulphides();
            break;

            case InstructionTypeRecalculateFFT:
            options->recalculateFFT();
            break;

            case InstructionTypeSetOutputDir:
            FileReader::setOutputDirectory(_outputDir);
            setMessage("Set output directory " + _outputDir);
            break;

            case InstructionTypeSetObjectValue:
            Notifiable::performObjectSet();
            break;
 
            case InstructionTypeGetObjectValue:
            Notifiable::performObjectGet();
            break;

            case InstructionTypePreviousState:
			options->previousState();
            break;

            default:
            break;
        }
        
		InstructionThread *thread = NULL;
		thread =  static_cast<InstructionThread *>(QThread::currentThread());

		if (thread->shouldDie())
		{
			QThread::currentThread()->exit(0);
		}

		display->setFocus();
        mutex.unlock();
    }
}

bool VagWindow::isRunningSomething()
{
    bool locked = mutex.try_lock();

    if (locked)
    {
        mutex.unlock();
        return false;
    }
    else
    {
        return true;
    }
}

void VagWindow::getPolymerMonomerCrystal(PolymerPtr *poly, CrystalPtr *cryst, MonomerPtr *monomer)
{
    OptionsPtr options = Options::getRuntimeOptions();
    *cryst = options->getActiveCrystal();

    *monomer = *(static_cast<MonomerPtr *>(getObject())); 
	*poly = (*monomer)->getPolymer();
}

void VagWindow::sidechainsToEnd()
{
	PolymerPtr polymer;	CrystalPtr crystal; MonomerPtr monomer;
	getPolymerMonomerCrystal(&polymer, &crystal, &monomer);
    polymer->refineToEnd(monomer->getResidueNum(), crystal, RefinementSidechain); 
}

void VagWindow::squeezeToEnd()
{
	PolymerPtr polymer;	CrystalPtr crystal; MonomerPtr monomer;
	getPolymerMonomerCrystal(&polymer, &crystal, &monomer);
    polymer->refineToEnd(monomer->getResidueNum(), crystal, RefinementModelRMSDZero); 
}

void VagWindow::modelPosToEnd()
{
	PolymerPtr polymer;	CrystalPtr crystal; MonomerPtr monomer;
	getPolymerMonomerCrystal(&polymer, &crystal, &monomer);
    polymer->refineToEnd(monomer->getResidueNum(), crystal, RefinementModelPos); 
}

void VagWindow::refineToEnd()
{
	PolymerPtr polymer;	CrystalPtr crystal; MonomerPtr monomer;
	getPolymerMonomerCrystal(&polymer, &crystal, &monomer);
    polymer->refineToEnd(monomer->getResidueNum(), crystal, RefinementFine); 
}

void VagWindow::enable()
{
	for (int i = 0; i < buttons.size(); i++)
	{
		buttons[i]->setEnabled(true);
	}
}

void VagWindow::disable()
{
	for (int i = 0; i < buttons.size(); i++)
	{
		buttons[i]->setEnabled(false);
	}
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

void VagWindow::pushFitWholeT()
{
    _instructionType = InstructionTypeFitWholeMoleculeTranslation;
    wait.wakeAll();
}

void VagWindow::pushFitWholeR()
{
    _instructionType = InstructionTypeFitWholeMoleculeRotation;
    wait.wakeAll();
}

void VagWindow::pushRefineFlexibility()
{
    _instructionType = InstructionTypeRefineFlexibility;
    wait.wakeAll();
}

void VagWindow::pushRefineDensity()
{
    _instructionType = InstructionTypeRefineDensity;
    wait.wakeAll();
}

void VagWindow::recalculateFFT()
{
    _instructionType = InstructionTypeRecalculateFFT;
    wait.wakeAll();
}

void VagWindow::pushBackboneAnalysis()
{
    _instructionType = InstructionTypeBackboneAnalysis;
    wait.wakeAll();
}

void VagWindow::findDisulphides()
{
    _instructionType = InstructionTypeFindDisulphides;
    wait.wakeAll();
}

void VagWindow::restorePreviousState()
{
    _instructionType = InstructionTypePreviousState;
	wait.wakeAll();
}

void VagWindow::pushExploreMcule()
{
    delete _explorer;
    _explorer = NULL;

    QObject *sent = sender();
    ChainMenuAction *action = static_cast<ChainMenuAction *>(sent);
    
    OptionsPtr options = Options::getRuntimeOptions();
    int crystalCount = options->crystalCount();
    
    MoleculePtr molecule = action->getMolecule();

    if (molecule)
    {
        _explorer = new MoleculeExplorer(this, molecule);
        _explorer->setGLKeeper(display->getKeeper());
        _explorer->show();
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

void VagWindow::openModel()
{
    delete _fileDialogue;
    _fileDialogue = new QFileDialog(this, tr("Open model"), tr("Model files (*.pdb, *.vbond)"));
    _fileDialogue->setFileMode(QFileDialog::AnyFile);
    _fileDialogue->show();
    
    QStringList fileNames;
    if (_fileDialogue->exec())
    {
        fileNames = _fileDialogue->selectedFiles();
    }
    
    if (fileNames.size() >= 1)
    {
        _pdbName = fileNames[0].toStdString();
        _instructionType = InstructionTypeOpenPDB;
        wait.wakeAll();
    }
}

void VagWindow::openMTZ()
{
    delete _fileDialogue;
    _fileDialogue = new QFileDialog(this, tr("Open MTZ"), tr("MTZ Files (*.mtz)"));
    _fileDialogue->setFileMode(QFileDialog::AnyFile);
    _fileDialogue->show();
    
    QStringList fileNames;
    if (_fileDialogue->exec())
    {
        fileNames = _fileDialogue->selectedFiles();
    }
    
    std::cout << "Understood " << fileNames.size() << std::endl;
    
    if (fileNames.size() >= 1)
    {
        _mtzName = fileNames[0].toStdString();
        _instructionType = InstructionTypeOpenMTZ;
        wait.wakeAll();
    }
    
}

void VagWindow::setOutput()
{
    delete _fileDialogue;
    _fileDialogue = new QFileDialog(this, tr("Set output directory"), tr("Folders"));
    _fileDialogue->setFileMode(QFileDialog::DirectoryOnly);
    _fileDialogue->show();
    
    QStringList fileNames;
    if (_fileDialogue->exec())
    {
        fileNames = _fileDialogue->selectedFiles();
    }
    
    std::cout << "Understood " << fileNames.size() << std::endl;
    
    if (fileNames.size() >= 1)
    {
        _outputDir = fileNames[0].toStdString();
        _instructionType = InstructionTypeSetOutputDir;
        wait.wakeAll();
    }
}

void VagWindow::wakeup()
{
    wait.wakeAll();
}

void VagWindow::setMessage(std::string message)
{
    Notifiable::setMessage(message);
    _lStatus->setText(QString::fromStdString(message));
}

VagWindow::~VagWindow()
{
    delete bSuperimpose;
    delete bRefinePos;
//    delete bRefineFlex;
    delete bExploreMolecule;
    delete bChangeBMult;
    delete bRecalculate;
    delete display;
    delete _myDialogue;
    delete _fileDialogue;
    delete _explorer;
}
