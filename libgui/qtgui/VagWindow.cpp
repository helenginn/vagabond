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

#ifdef __APPLE__
#define MENU_HAS_HEIGHT 0
#else
#define MENU_HAS_HEIGHT 1
#endif

void VagWindow::makeMenu()
{
	QMenu *scaling = menuBar()->addMenu(tr("&Scaling"));
	menus.push_back(scaling);
	
	_qaShell = scaling->addAction(tr("Shell-by-shell scaling"));
	_qaShell->setCheckable(true);
	connect(_qaShell, &QAction::triggered,
			[=]{ toggleScaling(ScalingTypeShell); });
	
	_qaKB = scaling->addAction(tr("Absolute scale + B factor"));
	_qaKB->setCheckable(true);
	connect(_qaKB, &QAction::triggered,
			[=]{ toggleScaling(ScalingTypeAbsBFactor); });
	
	_qaK = scaling->addAction(tr("Absolute scale only"));
	_qaK->setCheckable(true);
	connect(_qaK, &QAction::triggered,
			[=]{ toggleScaling(ScalingTypeAbs); });
	
	actions.push_back(_qaShell);
	actions.push_back(_qaKB);
	actions.push_back(_qaK);
	
	displayScaling();
	
	QAction *sep = scaling->addSeparator();
	actions.push_back(sep);
	
	QAction *adjust = scaling->addAction(tr("Adjust real-space B factor"));
	connect(adjust, &QAction::triggered, this, &VagWindow::adjustBFactor);

	QAction *chooseB = scaling->addAction(tr("Change real-space B factor..."));
	connect(chooseB, &QAction::triggered,
			[=]{ dialogueModify(Options::setGlobalBFactor, 
			                    "B factor applied in real space"); });
	actions.push_back(adjust);

	QMenu *model = menuBar()->addMenu(tr("&Model"));
	menus.push_back(model);
	
	QAction *samples = model->addAction(tr("Change model sampling..."));
	connect(samples, &QAction::triggered,
			[=]{ dialogueModify(Options::setNSamples, 
			                    "No. samples in ensemble:"); });
	
	actions.push_back(samples);
}

void VagWindow::dialogueModify(SimpleSet set, std::string title)
{
	delete _myDialogue;
	_myDialogue = new Dialogue(this, "New value", title, "100", "Set value");
	_myDialogue->setFunction(set);
	_myDialogue->setWindow(this);
	_myDialogue->show();

}

void VagWindow::makeButtons()
{
	buttons.clear();

	bRefinePos = new QPushButton("Refine positions to PDB", this);
	bRefinePos->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 0, BUTTON_WIDTH , 50);
	bRefinePos->setEnabled(false);
	connect(bRefinePos, &QPushButton::clicked,
			[=]{ pushSendInstruction(InstructionTypeRefinePositions); });
	buttons.push_back(bRefinePos);

	bFitWholeT = new QPushButton("Fit molecule translation", this);
	bFitWholeT->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 100, BUTTON_WIDTH , 50);
	bFitWholeT->setEnabled(false);
	connect(bFitWholeT, &QPushButton::clicked,
			[=]{ pushSendInstruction(InstructionTypeFitTranslation); });

	buttons.push_back(bFitWholeT);

	bRefineDensity = new QPushButton("Refine sidechains to density", this);
	bRefineDensity->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 200, BUTTON_WIDTH , 50);
	bRefineDensity->setEnabled(false);
	connect(bRefineDensity, &QPushButton::clicked,
			[=]{ pushSendInstruction(InstructionTypeRefineDensity); });
	buttons.push_back(bRefineDensity);

	bChangeBMult = new QPushButton("Set hetatm B multiplier", this);
	bChangeBMult->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 250, BUTTON_WIDTH , 50);
	bChangeBMult->setEnabled(false);
	connect(bChangeBMult, SIGNAL(clicked()), this, SLOT(pushBMultiplier()));
	buttons.push_back(bChangeBMult);

	bRecalculate = new QPushButton("Recalculate FFT", this);
	bRecalculate->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 300, BUTTON_WIDTH , 50);
	bRecalculate->setEnabled(false);
	connect(bRecalculate, &QPushButton::clicked,
			[=]{ pushSendInstruction(InstructionTypeRecalculateFFT); });
	buttons.push_back(bRecalculate);

	bWaterNetwork = new QPushButton("Scan bond params", this);
	bWaterNetwork->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 350, BUTTON_WIDTH , 50);
	bWaterNetwork->setEnabled(false);
	connect(bWaterNetwork, SIGNAL(clicked()), this, SLOT(refineWaterNetwork()));
	buttons.push_back(bWaterNetwork);

	bExploreMolecule = new QPushButton("Explore molecule", this);
	bExploreMolecule->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 450, BUTTON_WIDTH , 50);
	bExploreMolecule->setEnabled(false);
	bExploreMolecule->setMenu(new QMenu(this));
	buttons.push_back(bExploreMolecule);

	bExploreCrystal = new QPushButton("Explore crystal", this);
	bExploreCrystal->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 500, 
	                             BUTTON_WIDTH , 50);
	bExploreCrystal->setEnabled(false);
	connect(bExploreCrystal, SIGNAL(clicked()), this, 
	        SLOT(pushExploreCrystal()));
	buttons.push_back(bExploreCrystal);

	/*
	bPrevious = new QPushButton("Restore previous state", this);
	bPrevious->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 450, BUTTON_WIDTH , 50);
	bPrevious->setEnabled(false);
	connect(bPrevious, SIGNAL(clicked()), this, SLOT(restorePreviousState()));
	buttons.push_back(bPrevious);
	*/

	bCoot = new QPushButton("Open in Coot", this);
	bCoot->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 550, BUTTON_WIDTH , 50);
	bCoot->setEnabled(false);
	connect(bCoot, &QPushButton::clicked,
			[=]{ pushSendInstruction(InstructionTypeOpenInCoot); });
	buttons.push_back(bCoot);

	_myDialogue = NULL;
	_moleExplorer = NULL;
	_xtalExplorer = NULL;
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

	double menuHeight = 0;
	
	if (MENU_HAS_HEIGHT)
	{
		menuHeight += 30;
	}
	
	display->setGeometry(0, menuHeight, w - BUTTON_WIDTH, 
	                     h - STATUS_HEIGHT - menuHeight);
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
	display->setGeometry(0, 0, DEFAULT_WIDTH - BUTTON_WIDTH, 
	                     DEFAULT_HEIGHT - STATUS_HEIGHT);

	makeMenu();
	makeButtons();

	QFont bigFont = QFont("Helvetica", 16);
	_lStatus = new QLabel("Vagabond at your service.", this);
	_lStatus->setFont(bigFont);
	_lStatus->setGeometry(10, DEFAULT_HEIGHT - STATUS_HEIGHT, 680, STATUS_HEIGHT);
	_lStatus->show();

	_argc = argc;
	_argv = argv;
	display->setFocus();
	display->setFocusPolicy(Qt::StrongFocus);

	_instructionType = InstructionTypeNone;
	_instructionThread.setVagWindow(this);
	_instructionThread.start(); 
}

int VagWindow::waitForInstructions()
{
	OptionsPtr options = Options::getRuntimeOptions();
	options->setNotify((Notifiable *)this);
	InstructionThread *thread = NULL;
	thread =  static_cast<InstructionThread *>(QThread::currentThread());

	try
	{
		options->run();
	}
	catch (int e)
	{
		if (e == 10)
		{
			return 1;
		}
	}

	while (true)
	{
		try
		{
			/* Returned instructions to GUI */
			switch (_instructionType)
			{
				case InstructionTypeResetExplorer:
				updateExplorerButton();
				break;

				default:
				break;
			}

			mutex.lock();
			wait.wait(&mutex);

			/* GUI instructions to Vagabond */
			switch (_instructionType)
			{
				case InstructionTypeRefinePositions:
				options->refineAll(RefinementModelPos, 1);
				break;

				case InstructionTypeFitTranslation: 
				options->fitWholeMolecule(true, false);
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

				case InstructionTypeFindDisulphides:
				options->findDisulphides();
				break;

				case InstructionTypeRefineWaterNetwork:
				options->scanBondParams();
//				options->refineAll(RefinementWaterNetwork, 1);
				break;

				case InstructionTypeRecalculateFFT:
				options->recalculateFFT();
				break;

				case InstructionTypeOpenInCoot:
				options->openInCoot();
				break;

				case InstructionTypeSetObjectValue:
				Notifiable::performObjectSet();
				break;

				case InstructionTypeGetObjectValue:
				Notifiable::performObjectGet();
				_moleExplorer->updateCorrelation();
				break;

				case InstructionTypePreviousState:
				options->previousState();
				break;
				
				case InstructionTypeSplitBond:
				splitBond();
				break;
				
				case InstructionTypeAdjustBFactor:
				options->adjustBFactor();
				break;

				default:
				break;
			}
		}
		catch (int e)
		{
			if (e == 10)
			{
				return 1;
			}
		}

		if (thread->shouldDie())
		{
			return 0;
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

void VagWindow::pushSendInstruction(InstructionType inst)
{
	_instructionType = inst;
	wait.wakeAll();
}

void VagWindow::refineWaterNetwork()
{
	_instructionType = InstructionTypeRefineWaterNetwork;
	wait.wakeAll();
}

void VagWindow::restorePreviousState()
{
	_instructionType = InstructionTypePreviousState;
	wait.wakeAll();
}

void VagWindow::pushExploreCrystal()
{
	delete _xtalExplorer;
	_xtalExplorer = NULL;
	
	OptionsPtr options = Options::getRuntimeOptions();
	CrystalPtr crystal = options->getActiveCrystal();
	
	if (crystal)
	{
		_xtalExplorer = new CrystalExplorer(this, crystal);
		_xtalExplorer->show();
	}
}

void VagWindow::pushExploreMcule()
{
	QObject *sent = sender();
	ChainMenuAction *action = static_cast<ChainMenuAction *>(sent);

	OptionsPtr options = Options::getRuntimeOptions();
	int crystalCount = options->crystalCount();

	MoleculePtr molecule = action->getMolecule();

	if (molecule)
	{
		delete _moleExplorer;
		_moleExplorer = NULL;

		_moleExplorer = new MoleculeExplorer(this, molecule);
		_moleExplorer->setGLKeeper(display->getKeeper());
		_moleExplorer->show();
	}
}

void VagWindow::splitBond()
{
	double num = getValue();
	void *obj = getObject();
	Bond *bond = static_cast<Bond *>(obj);

	Bond *blockade = bond;

	/* Find split blockade */
	for (int i = 0; i < num; i++)
	{
		if (blockade->downstreamBondGroupCount() && 
		    blockade->downstreamBondCount(0))
		{
			AtomPtr newAtom = blockade->downstreamAtom(0, 0);
			blockade = &*(ToBondPtr(newAtom->getModel())); 
		}
	}
	
	if (num >= 0)
	{
		std::cout << "Blocking split at " << blockade->shortDesc() 
		<< std::endl;
		blockade->setSplitBlock();
	}
	
	std::cout << "Splitting bond at " << bond->shortDesc() << std::endl;
	bond->splitBond();

	bond->getMinor()->getMolecule()->refreshPositions();
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

void VagWindow::wakeup()
{
	wait.wakeAll();
}

void VagWindow::displayScaling()
{
	ScalingType type = Options::getScalingType();

	_qaShell->setChecked(type == ScalingTypeShell);
	_qaKB->setChecked(type == ScalingTypeAbsBFactor);
	_qaK->setChecked(type == ScalingTypeAbs);
}

void VagWindow::toggleScaling(ScalingType type)
{
	Options::setScalingType(type);
	displayScaling();
}

void VagWindow::setMessage(std::string message)
{
	Notifiable::setMessage(message);
	_lStatus->setText(QString::fromStdString(message));
}

void VagWindow::setRenderDensity()
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	display->renderDensity(crystal);
}

void VagWindow::adjustBFactor()
{
	_instructionType = InstructionTypeAdjustBFactor;
	wait.wakeAll();
}

VagWindow::~VagWindow()
{
	delete bRefinePos;
	delete bExploreMolecule;
	delete bChangeBMult;
	delete bRecalculate;
	delete display;
	delete _myDialogue;
	delete _fileDialogue;
	delete _moleExplorer;
}


