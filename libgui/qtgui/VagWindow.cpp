//
//  VagWindow.cpp
//  Vagabond
//
//  Created by Helen Ginn on 11/12/2016.
//  Copyright (c) 2016 Helen Ginn. All rights reserved.
//

#define DEFAULT_WIDTH 900
#define DEFAULT_HEIGHT 660

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
#include "ChainMenuAction.h"
#include "../../libsrc/FileReader.h"

void VagWindow::makeMenu()
{
    _fileDialogue = NULL;

    QMenu *fileMenu = menuBar()->addMenu(tr("&File"));
    QAction *aOutput = fileMenu->addAction(tr("Set &output"));
    connect(aOutput, &QAction::triggered, this, &VagWindow::setOutput);
    QAction *aOpenPDB = fileMenu->addAction(tr("Open &PDB..."));
    connect(aOpenPDB, &QAction::triggered, this, &VagWindow::openPDB);
    QAction *aOpenMTZ = fileMenu->addAction(tr("Open &MTZ..."));
    connect(aOpenMTZ, &QAction::triggered, this, &VagWindow::openMTZ);
}

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

    bRefineDensity = new QPushButton("Refine sidechains to density", this);
    bRefineDensity->setGeometry(700, 150, 200, 50);
    bRefineDensity->setEnabled(false);
    connect(bRefineDensity, SIGNAL(clicked()), this, SLOT(pushRefineDensity()));

    bRecalculate = new QPushButton("Recalculate FFT", this);
    bRecalculate->setGeometry(700, 200, 200, 50);
    bRecalculate->setEnabled(false);
    connect(bRecalculate, SIGNAL(clicked()), this, SLOT(recalculateFFT()));    
    bChangeBMult = new QPushButton("Set hetatm B multiplier", this);
    bChangeBMult->setGeometry(700, 250, 200, 50);
    bChangeBMult->setEnabled(false);
    connect(bChangeBMult, SIGNAL(clicked()), this, SLOT(pushBMultiplier()));
    
    bExploreMolecule = new QPushButton("Explore molecule", this);
    bExploreMolecule->setGeometry(700, 300, 200, 50);
    bExploreMolecule->setEnabled(false);
    bExploreMolecule->setMenu(new QMenu(this));

    
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

    CrystalPtr crystal = options->getCrystal(0);

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

VagWindow::VagWindow(QWidget *parent,
                     int argc, char *argv[]) : QMainWindow(parent)
{
    this->resize(DEFAULT_WIDTH, DEFAULT_HEIGHT);
    
    this->setWindowTitle("Vagabond");
    
    display = new VagabondGLWidget(this);
    display->setGeometry(0, 0, 700, 600);
    
    makeButtons();
    makeMenu();    

    QFont bigFont = QFont("Helvetica", 16);
    _lStatus = new QLabel("Vagabond at your service.", this);
    _lStatus->setFont(bigFont);
    _lStatus->setGeometry(10, 600, 680, 60);
    _lStatus->show();

    _argc = argc;
    _argv = argv;
    
    _instructionType = InstructionTypeNone;
    _instructionThread.setVagWindow(this);
    _instructionThread.start(); 
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

            case InstructionTypeRefineDensity: 
            options->refineAll(RefinementFine, 1);
            break;

            case InstructionTypeChangeBMult:
            options->applyBMultiplier();
            break;
       
            case InstructionTypeOpenPDB:
            options->openPDB(_pdbName);
            updateExplorerButton();
            break;

            case InstructionTypeOpenMTZ:
            options->openMTZ(_mtzName);
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

            default:
            break;
        }
        
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

void VagWindow::enable()
{
    bSuperimpose->setEnabled(true);
    bRefinePos->setEnabled(true);
    bRefineFlex->setEnabled(true);
    bChangeBMult->setEnabled(true);
    bExploreMolecule->setEnabled(true);
    bRecalculate->setEnabled(true);
    bRefineDensity->setEnabled(true);
}

void VagWindow::disable()
{
    bSuperimpose->setEnabled(false);
    bRefinePos->setEnabled(false);
    bRefineFlex->setEnabled(false);
    bChangeBMult->setEnabled(false);
    bExploreMolecule->setEnabled(false);
    bRecalculate->setEnabled(false);
    bRefineDensity->setEnabled(false);
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

void VagWindow::openPDB()
{
    delete _fileDialogue;
    _fileDialogue = new QFileDialog(this, tr("Open PDB"), tr("PDB Files (*.pdb)"));
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
    delete bRefineFlex;
    delete bExploreMolecule;
    delete bChangeBMult;
    delete bRecalculate;
    delete display;
    delete _myDialogue;
    delete _fileDialogue;
    delete _explorer;
}
