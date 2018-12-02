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

void VagWindow::menuItem(QMenu *menu, std::string title,
                         InstructionType instr)
{
	QAction *action = menu->addAction(QString::fromStdString(title));
	connect(action, &QAction::triggered,
			[=]{ pushSendInstruction(instr); });
	actions.push_back(action);
}

void VagWindow::makeMenu()
{
	QMenu *scaling = menuBar()->addMenu(tr("&Scaling"));
	menus.push_back(scaling);
	
	_qaShell = scaling->addAction(tr("Shell-by-shell scaling"));
	_qaShell->setCheckable(true);
	connect(_qaShell, &QAction::triggered,
			[=]{ toggleScaling(ScalingTypeShell); });
	
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
	
	QAction *bmult = model->addAction(tr("B factor HETATM multiplier..."));
	connect(bmult, &QAction::triggered,
			[=]{ dialogueModify(Options::setBMult, 
			                    "Multiply original B factor of HETATMs by...",
		                     1.0); });
	actions.push_back(bmult);
	
	QAction *bsubt = model->addAction(tr("B factor HETATM subtract..."));
	connect(bsubt, &QAction::triggered,
			[=]{ dialogueModify(Options::setBSubt, 
			                    "Subtract from B factor of HETATMs by...",
		                     0.0); });
	actions.push_back(bsubt);
	
	menuItem(model, "Omit scan", InstructionTypeOmitScan);
	menuItem(model, "Find flexibility", InstructionTypeReflex);
	
	QAction *refit = model->addAction(tr("Refit backbone region..."));
	connect(refit, SIGNAL(triggered()), this, 
	        SLOT(refitBackbone()));
	actions.push_back(refit);

	sep = model->addSeparator();
	actions.push_back(sep);

	QAction *expMol = model->addAction(tr("Explore crystal..."));
	connect(expMol, SIGNAL(triggered()), this, 
	        SLOT(pushExploreCrystal()));
	actions.push_back(bsubt);

	menuItem(model, "Open in Coot", InstructionTypeOpenInCoot);

	QMenu *mRefine = menuBar()->addMenu(tr("&Refine"));
	menus.push_back(mRefine);

	menuItem(mRefine, "Positions to PDB", 
	         InstructionTypeRefinePositions);
	menuItem(mRefine, "Intermolecule movements",
	         InstructionTypeFitTranslation);
	menuItem(mRefine, "Intramolecule movements",
	         InstructionTypeRefineIntramolecule);
	menuItem(mRefine, "Sidechains to density",
	         InstructionTypeRefineDensity);
	menuItem(mRefine, "Recalculate FFT",
	         InstructionTypeRecalculateFFT);
}

void VagWindow::refitBackbone()
{
	delete _myDialogue;
	_myDialogue = new Dialogue(this, "Specify backbone region", 
	                           "Two limiting residues", 
	                           "795 805",
	                           "Refit");
	_myDialogue->setWindow(this);
	_myDialogue->setTag(DialogueRefit);
	_myDialogue->show();
}

void VagWindow::dialogueModify(SimpleSet set, std::string title,
                               double _default)
{
	delete _myDialogue;
	_myDialogue = new Dialogue(this, "New value", title, 
	                           f_to_str(_default, 0),
	                           "Set value");
	_myDialogue->setFunction(set);
	_myDialogue->setWindow(this);
	_myDialogue->show();

}

void VagWindow::makeButtons()
{
	buttons.clear();

	/*
	bChelate = new QPushButton("Chelate metals", this);
	bChelate->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 400, BUTTON_WIDTH , 50);
	bChelate->setEnabled(false);
	connect(bChelate, &QPushButton::clicked,
			[=]{ pushSendInstruction(InstructionTypeChelate); });
	buttons.push_back(bChelate);
	*/

	bExploreMolecule = new QPushButton("Explore molecule", this);
	bExploreMolecule->setGeometry(DEFAULT_WIDTH - BUTTON_WIDTH, 450, BUTTON_WIDTH , 50);
	bExploreMolecule->setEnabled(false);
	bExploreMolecule->setMenu(new QMenu(this));
	buttons.push_back(bExploreMolecule);

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

	CrystalPtr crystal = options->getActiveCrystal();

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
				crystal->refinePositions();
				options->recalculateFFT();
				break;

				case InstructionTypeWhack: 
				options->whack();
				break;

				case InstructionTypeFitTranslation: 
				crystal->fitWholeMolecules();
				options->recalculateFFT();
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
				crystal->refineSidechains();
				options->recalculateFFT();
				break;

				case InstructionTypeChangeBMult:
				options->applyBMultiplier();
				break;

				case InstructionTypeFindDisulphides:
				options->findDisulphides();
				break;

				case InstructionTypeRefineIntramolecule:
				crystal->refineIntraMovements();
				options->recalculateFFT();
				break;

				case InstructionTypeRecalculateFFT:
				options->recalculateFFT();
				break;

				case InstructionTypeChelate:
				options->chelate();
				break;

				case InstructionTypeOpenInCoot:
				options->openInCoot();
				break;

				case InstructionTypeReflex:
				options->reflex();
				break;

				case InstructionTypeOmitScan:
				options->omitScan();
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
				
				case InstructionTypeRefitBackbone:
				options->refitBackbone(_rangeStart, _rangeEnd);
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
		_xtalExplorer->setKeeper(display->getKeeper());
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

void VagWindow::pause(bool on)
{
	display->getKeeper()->pause(on);
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

	if (type == DialogueRefit)
	{
		_rangeStart = trial[0];
		_rangeEnd = trial[1];
		_instructionType = InstructionTypeRefitBackbone;
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
//	_qaKB->setChecked(type == ScalingTypeAbsBFactor);
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

void VagWindow::focusOnPosition(vec3 pos)
{
	display->getKeeper()->focusOnPosition(pos);
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


