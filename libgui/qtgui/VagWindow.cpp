//
//  VagWindow.cpp
//  Vagabond
//
//  Created by Helen Ginn on 11/12/2016.
//  Copyright (c) 2016 Helen Ginn. All rights reserved.
//

#define DEFAULT_WIDTH 900
#define DEFAULT_HEIGHT 660
#define LOG_WIDTH 350
#define STATUS_HEIGHT 60
#define BUTTON_WIDTH 0

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
#include "../../libsrc/charmanip.h"
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
	actions.push_back(_qaK);
	
	displayScaling();
	
	QAction *sep = scaling->addSeparator();
	actions.push_back(sep);
	
	/*
	QAction *adjust = scaling->addAction(tr("Adjust real-space B factor"));
	connect(adjust, &QAction::triggered, this, &VagWindow::adjustBFactor);
	actions.push_back(adjust);
	*/

	QAction *chooseB = scaling->addAction(tr("Change real-space B factor..."));
	connect(chooseB, &QAction::triggered,
			[=]{ dialogueModify(Options::setGlobalBFactor, 
			                    "B factor applied in real space"); });
	actions.push_back(chooseB);

	QMenu *model = menuBar()->addMenu(tr("&Model"));
	menus.push_back(model);
	
	QAction *samples = model->addAction(tr("Change model sampling..."));
	connect(samples, &QAction::triggered,
			[=]{ dialogueModify(Options::changeSamplesAndFit, 
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
	
	QAction *bprobe = model->addAction(tr("Solvent probe radius..."));
	connect(bprobe, &QAction::triggered,
			[=]{ dialogueModify(Options::setProbeRadius, 
			                    "Set solvent probe radius (Å)",
		                     0.0); });
	actions.push_back(bprobe);
	
	/*
	menuItem(model, "Omit scan", InstructionTypeOmitScan);
	menuItem(model, "Find flexibility", InstructionTypeReflex);
	
	QAction *refit = model->addAction(tr("Refit backbone region..."));
	connect(refit, SIGNAL(triggered()), this, 
	        SLOT(refitBackbone()));
	actions.push_back(refit);
	
	refit = model->addAction(tr("Fix erroneous zones..."));
	connect(refit, SIGNAL(triggered()), this, 
	        SLOT(fixErroneousZones()));
	actions.push_back(refit);
	*/

	sep = model->addSeparator();
	actions.push_back(sep);

	QAction *expMol = model->addAction(tr("Explore crystal..."));
	connect(expMol, SIGNAL(triggered()), this, 
	        SLOT(pushExploreCrystal()));
	actions.push_back(bsubt);

	QAction *coot = model->addAction("Open in Coot");
	connect(coot, SIGNAL(triggered()), this, 
	        SLOT(openInCoot()));
//	actions.push_back(coot);

	QMenu *mRefine = menuBar()->addMenu(tr("&Refine"));
	menus.push_back(mRefine);

	menuItem(mRefine, "Positions to PDB", 
	         InstructionTypeRefinePositions);
	menuItem(mRefine, "Backbone to density", 
	         InstructionTypeRefinePosToDensity);
	menuItem(mRefine, "Intermolecule movements",
	         InstructionTypeFitTranslation);
	menuItem(mRefine, "Intramolecule movements",
	         InstructionTypeRefineIntramolecule);
	menuItem(mRefine, "Sidechain positions to density", 
	         InstructionTypeRefineSidePos);
	menuItem(mRefine, "Sidechains to density",
	         InstructionTypeRefineDensity);
	menuItem(mRefine, "Recalculate FFT",
	         InstructionTypeRecalculateFFT);

	QAction *action = mRefine->addAction("Cancel further refinement");
	connect(action, &QAction::triggered,
			[=]{ setShouldCancel(); });
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

void VagWindow::gotoResidueDialogue()
{
	delete _myDialogue;
	_myDialogue = new Dialogue(this, "Go to residue", 
	                           "Choose chain and number",
	                           "A320",
	                           "Go");
	_myDialogue->setWindow(this);
	_myDialogue->setTag(DialogueGoto);
	_myDialogue->show();
}

void VagWindow::dialogueModify(Setter set, std::string title,
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

	_myDialogue = NULL;
	_xtalExplorer = NULL;
	_errorExplorer = NULL;
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
	
	double logWidth = LOG_WIDTH;
	if (w * 0.3 > logWidth)
	{
		logWidth = w * 0.3;
	}
	double displWidth = (!_showingLog ? w : w - logWidth);

	display->setGeometry(0, menuHeight, displWidth, 
	                     h - STATUS_HEIGHT - menuHeight);
	display->resizeGL();

	_lStatus->setGeometry(10, h - STATUS_HEIGHT, w - BUTTON_WIDTH, 
	                      STATUS_HEIGHT);
	
	if (_showingLog)
	{
		_logView->setGeometry(w - logWidth, menuHeight - 1, logWidth, 
		                      h - STATUS_HEIGHT - menuHeight);
		_logView->show();
	}
	else
	{
		_logView->hide();
	}

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
	display->setVagWindow(this);
	display->setGeometry(0, 0, DEFAULT_WIDTH - BUTTON_WIDTH, 
	                     DEFAULT_HEIGHT - STATUS_HEIGHT);

	makeMenu();
	makeButtons();

	QFont bigFont = QFont("Helvetica", 16);
	_lStatus = new QLabel("Vagabond at your service.", this);
	_lStatus->setFont(bigFont);
	_lStatus->setGeometry(10, DEFAULT_HEIGHT - STATUS_HEIGHT, 680, STATUS_HEIGHT);
	_lStatus->show();
	
	QFont monoFont("Monospace");
	monoFont.setPointSize(8);
	monoFont.setStyleHint(QFont::Monospace);
	_logView = new QTextEdit(this);
	_logView->setAcceptRichText(false);
	_logView->setFont(monoFont);
	_logView->setReadOnly(true);
	_logView->setGeometry(DEFAULT_WIDTH - LOG_WIDTH, 0, LOG_WIDTH, 
	                      DEFAULT_WIDTH - STATUS_HEIGHT);

	_argc = argc;
	_argv = argv;
	_showingLog = false;
	display->setFocus();
	display->setFocusPolicy(Qt::StrongFocus);


	_instructionType = InstructionTypeNone;
	_instructionThread.setVagWindow(this);

	connect(this, SIGNAL(errorReceived(Shouter *)), 
	        this, SLOT(displayMessage(Shouter *)), Qt::QueuedConnection);
	
	connect(this, SIGNAL(appendSignal()),
	        this, SLOT(append()), Qt::UniqueConnection);

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
	catch (Shouter *shout)
	{
		emit errorReceived(shout);
		options->clear();
		return 1;
	}
	
	_startScreen->finishUp();

	CrystalPtr crystal = options->getActiveCrystal();

	while (true)
	{
		try
		{
			/* Returned instructions to GUI */
			switch (_instructionType)
			{
				case InstructionTypeResetExplorer:
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
				disable();
				crystal->refinePositions();
				options->recalculateFFT();
				enable();
				break;
				
				case InstructionTypeRefinePosToDensity:
				disable();
				crystal->refineCrude();
				options->recalculateFFT();
				enable();
				break;

				case InstructionTypeFitTranslation: 
				crystal->fitWholeMolecules();
				options->recalculateFFT();
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
				disable();
				crystal->refineSidechains();
				options->recalculateFFT();
				enable();
				break;

				case InstructionTypeFindDisulphides:
				options->findDisulphides();
				break;
				
				case InstructionTypeRefineSidePos:
				disable();
				crystal->refineSidechainPositions();
				options->recalculateFFT();
				enable();
				break;

				case InstructionTypeRefineIntramolecule:
				disable();
				crystal->refineIntraMovements();
				options->recalculateFFT();
				enable();
				break;

				case InstructionTypeRecalculateFFT:
				options->recalculateFFT();
				break;

				case InstructionTypeChelate:
				options->chelate();
				break;

				case InstructionTypeOmitScan:
				options->omitScan();
				break;

				case InstructionTypeSetObjectValue:
				Notifiable::performObjectSet();
				break;

				case InstructionTypeGetObjectValue:
				Notifiable::performObjectGet();
				_xtalExplorer->updateCorrelation();
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

				case InstructionTypeManualRefine:
				display->manualRefine();
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
	if (!Notifiable::isEnabled())
	{
		return true;
	}

	bool locked = mutex.tryLock();

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
	for (int i = 0; i < actions.size(); i++)
	{
		actions[i]->setEnabled(true);
	}
	
	Notifiable::enable();
}

void VagWindow::disable()
{
	for (int i = 0; i < actions.size(); i++)
	{
		actions[i]->setEnabled(false);
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
		_xtalExplorer->setVagWindow(this);
		_xtalExplorer->setKeeper(display->getKeeper());
		_xtalExplorer->show();
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

void VagWindow::receiveGotoResidue(std::string diagString)
{
	if (diagString.length() == 0)
	{
		std::cout << "No residue entered." << std::endl;
		return;
	}

	to_upper(diagString);

	char *start = &diagString[0];
	std::string chain = "";

	int resNum = 0;
	
	if (*start < '0' || *start > '9')
	{
		chain += *start;
		
		if (diagString.length() <= 1)
		{
			resNum = -INT_MAX;
		}
		else
		{
			start++;
			resNum = atoi(start);
		}
	}
	
	display->getKeeper()->selectResidue(chain, resNum);

	_myDialogue->hide();
	_myDialogue->deleteLater();
	_myDialogue = NULL;
}

void VagWindow::receiveDialogue(DialogueType type, std::string diagString)
{
	if (type == DialogueGoto)
	{
		receiveGotoResidue(diagString);
		return;
	}

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
	_myDialogue->deleteLater();
	_myDialogue = NULL;
}

void VagWindow::wakeup()
{
	wait.wakeAll();
}

void VagWindow::displayMessage(Shouter *shout)
{
	QString message = QString::fromStdString(shout->getMessage());
	QMessageBox::critical(this, "Problem", message, QMessageBox::Ok);
	this->hide();
	delete shout;
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

void VagWindow::fixErroneousZones()
{
	delete _errorExplorer;
	_errorExplorer = NULL;
	
	OptionsPtr options = Options::getRuntimeOptions();
	CrystalPtr crystal = options->getActiveCrystal();
	
	if (crystal)
	{
		_errorExplorer = new ErroneousZone(this, crystal);
		_errorExplorer->show();
	}
}

void VagWindow::openInCoot()
{
	Options::getActiveCrystal()->openInCoot();
}

void VagWindow::toggleLog()
{
	_showingLog = !_showingLog;
	resizeEvent(NULL);
}

void VagWindow::append()
{
	_logView->moveCursor(QTextCursor::End);
	_guiOutMut.lock();
	QString qOut = QString::fromStdString(_guiOut);
	_guiOut.clear();
	_guiOutMut.unlock();
	_logView->insertPlainText(qOut);
	_logView->moveCursor(QTextCursor::End);
}

void VagWindow::appendToLog(char *msg)
{
	if (*msg != NULL)
	{
		/* Will wait until main thread unlocks */
		_guiOutMut.lock();
		_guiOut += msg;
		_guiOutMut.unlock();

		emit appendSignal();
	}
}

VagWindow::~VagWindow()
{
	delete bRefinePos;
	delete bRecalculate;
	delete display;
	delete _myDialogue;
	delete _fileDialogue;
}


