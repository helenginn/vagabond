// Vagabond
// Copyright (C) 2017-2018 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#define DEFAULT_WIDTH (600)
#define DEFAULT_HEIGHT (450)
#define TEXTBOX_X (175)
#define LABEL_WIDTH (150)
#define CHECKBOX_WIDTH (250)
#define TEXTBOX_WIDTH (275)

#include "StartScreen.h"
#include "../../libsrc/Options.h"
#include "../../libsrc/Shouter.h"
#include "VagWindow.h"
#include <QtWidgets/qmessagebox.h>
#include "../../libsrc/FileReader.h"

std::string StartScreen::findNewFolder()
{
	int count = 1;
	
	while (1)
	{
		std::string test = "refine_" + i_to_str(count);

		if (!file_exists(test))
		{
			return test;
		}
		
		count++;
	}
}

void StartScreen::makeButtons()
{
	_bRun = new QPushButton("Start", this);
	_bRun->setGeometry(250, DEFAULT_HEIGHT - 50, 100, 25);
	_bRun->setFocus();
	connect(_bRun, SIGNAL(clicked()), this, SLOT(pushRun()));
	
	int height = 25;
	
	QLabel *q = new QLabel("Input files and output folders:", this);
	q->setGeometry(25, height, 300, 20);
	_widgets.push_back(q);
	 
	height += 40;

	q = new QLabel("Reflection file:", this);
	q->setGeometry(25, height, LABEL_WIDTH, 20);
	_widgets.push_back(q);

	std::string mtzFile = _options->_mtzFile;
	_tMtz = new QLineEdit(QString::fromStdString(mtzFile), this);
	_tMtz->setGeometry(TEXTBOX_X, height, TEXTBOX_WIDTH, 25);

	_bMtz = new QPushButton("Choose MTZ", this);
	_bMtz->setGeometry(475, height, 100, 25);
	connect(_bMtz, SIGNAL(clicked()), this, SLOT(chooseMtz()));

	height += 40;

	q = new QLabel("Model file:", this);
	q->setGeometry(25, height, LABEL_WIDTH, 20);
	_widgets.push_back(q);

	std::string modelFile = _options->_modelFile;
	_tModel = new QLineEdit(QString::fromStdString(modelFile), this);
	_tModel->setGeometry(TEXTBOX_X, height, TEXTBOX_WIDTH, 25);
	
	_bModel = new QPushButton("Choose model", this);
	_bModel->setGeometry(475, height, 100, 25);
	connect(_bModel, SIGNAL(clicked()), this, SLOT(chooseModel()));

	height += 40;

	q = new QLabel("Output folder:", this);
	q->setGeometry(25, height, LABEL_WIDTH, 20);
	_widgets.push_back(q);

	std::string outputDir = _options->_outputDir;
	
	if (outputDir.length() == 0)
	{
		outputDir = findNewFolder();
	}
	
	_tDir = new QLineEdit(QString::fromStdString(outputDir), this);
	_tDir->setGeometry(TEXTBOX_X, height, TEXTBOX_WIDTH, 25);
	
	_bDir = new QPushButton("Choose folder", this);
	_bDir->setGeometry(475, height, 100, 25);
	connect(_bDir, SIGNAL(clicked()), this, SLOT(chooseDir()));

	height += 50;
	int indent = 300;
	
	_showOpts = false;
	_bOptionals = new QPushButton("Show optional parameters", this);
	_bOptionals->setGeometry(indent, height, 250, 30);
	connect(_bOptionals, SIGNAL(clicked()), this, SLOT(toggleOptionals()));
	_widgets.push_back(_bOptionals);
	
	int top = height;
	height += 40;

	q = new QLabel("Min resolution (Å):", this);
	q->setGeometry(indent, height, LABEL_WIDTH, 20);
	q->hide();
	_widgets.push_back(q);
	_optWidgets.push_back(q);

	height += 40;

	q = new QLabel("Max resolution (Å):", this);
	q->setGeometry(indent, height, LABEL_WIDTH, 20);
	q->hide();
	_widgets.push_back(q);
	_optWidgets.push_back(q);
	
	indent += 150;
	height -= 40;

	_tMinRes = new QLineEdit("", this);
	_tMinRes->setGeometry(indent, height, 120, 25);
	_tMinRes->hide();
	_optWidgets.push_back(_tMinRes);
	
	if (Options::minRes() > 0)
	{
		double res = Options::minRes();
		std::string text = f_to_str(res, 3);
		_tMinRes->setText(QString::fromStdString(text));
	}

	height += 40;
	_tMaxRes = new QLineEdit("", this);
	_tMaxRes->setGeometry(indent, height, 120, 25);
	_tMaxRes->hide();
	_optWidgets.push_back(_tMaxRes);

	if (Options::maxRes() > 0)
	{
		double res = Options::maxRes();
		std::string text = f_to_str(res, 3);
		_tMaxRes->setText(QString::fromStdString(text));
	}

	height += 40;
	indent -= 150;

	q = new QLabel("Solvent add radius (Å):", this);
	q->setGeometry(indent, height, LABEL_WIDTH, 20);
	q->hide();
	_widgets.push_back(q);
	_optWidgets.push_back(q);

	indent += 150;
	
	_tRadius = new QLineEdit("", this);
	_tRadius->setGeometry(indent, height, 120, 25);
	_tRadius->hide();
	_optWidgets.push_back(_tRadius);

	if (Options::getProbeRadius() >= 0)
	{
		double radius = Options::getProbeRadius();
		std::string text = f_to_str(radius, 2);
		_tRadius->setText(QString::fromStdString(text));
	}

	height += 40;
	indent -= 150;

	q = new QLabel("Relative flex (AU):", this);
	q->setGeometry(indent, height, LABEL_WIDTH, 20);
	q->hide();
	_widgets.push_back(q);
	_optWidgets.push_back(q);

	indent += 150;
	
	_tRelFlex = new QLineEdit("", this);
	_tRelFlex->setGeometry(indent, height, 120, 25);
	_tRelFlex->hide();
	_optWidgets.push_back(_tRelFlex);

	if (Options::getBStart() >= 0)
	{
		double flex = Options::getBStart();
		std::string text = f_to_str(flex, 2);
		_tRelFlex->setText(QString::fromStdString(text));
	}

	/** Refinement options **/
	
	height = top;
	indent = 25;

	q = new QLabel("Refinement options", this);
	q->setGeometry(indent, height, CHECKBOX_WIDTH, 20);
	_widgets.push_back(q);
	
	height += 20;
	_cRefine = new QCheckBox("Perform refinement, including", this);
	_cRefine->setChecked(Options::_refine);
	_cRefine->setGeometry(indent, height, CHECKBOX_WIDTH, 20);
	connect(_cRefine, SIGNAL(stateChanged(int)), 
	        this, SLOT(toggleDisableOptions(int)));

	indent += 10;
	height += 20;
	_cPosition = new QCheckBox("positions to match PDB", this);
	_cPosition->setChecked(Options::_rPosition);
	_cPosition->setGeometry(indent, height, CHECKBOX_WIDTH, 20);
	_allToggle.push_back(_cPosition);

	height += 20;
	_cInter = new QCheckBox("whole molecule flex", this);
	_cInter->setChecked(Options::_rInter);
	_cInter->setGeometry(indent, height, CHECKBOX_WIDTH, 20);
	_allToggle.push_back(_cInter);

	height += 20;
	_cIntra = new QCheckBox("intramolecular flex", this);
	_cIntra->setChecked(Options::_rIntra);
	_cIntra->setGeometry(indent, height, CHECKBOX_WIDTH, 20);
	_allToggle.push_back(_cIntra);

	height += 20;
	_cCbAngles = new QCheckBox("variable Cb angles (from PDB)", this);
	_cCbAngles->setChecked(Options::_bondAngles > 0);
	_cCbAngles->setGeometry(indent, height, CHECKBOX_WIDTH, 20);
	_allToggle.push_back(_cCbAngles);
	connect(_cCbAngles, SIGNAL(stateChanged(int)), 
	        this, SLOT(toggleDisableOptions(int)));

	indent += 10;
	height += 20;
	_cCgAngles = new QCheckBox("... and Cg angles for aromatics", this);
	_cCgAngles->setChecked(Options::_bondAngles > 1);
	_cCgAngles->setGeometry(indent, height, CHECKBOX_WIDTH, 20);
	_allToggle.push_back(_cCgAngles);
	_angleToggle.push_back(_cCgAngles);
	connect(_cCgAngles, SIGNAL(stateChanged(int)), 
	        this, SLOT(toggleDisableOptions(int)));

	indent += 10;
	height += 20;
	_cGlyAngles = new QCheckBox("... and glycine backbone angles", this);
	_cGlyAngles->setChecked(Options::_bondAngles > 2);
	_cGlyAngles->setGeometry(indent, height, CHECKBOX_WIDTH, 20);
	_allToggle.push_back(_cGlyAngles);
	_angleToggle.push_back(_cGlyAngles);
	_angle2Toggle.push_back(_cGlyAngles);

	indent -= 20;
	height += 20;
	_cSidechains = new QCheckBox("refine sidechains to density", this);
	_cSidechains->setChecked(Options::_rSidechains);
	_cSidechains->setGeometry(indent, height, CHECKBOX_WIDTH, 20);
	_allToggle.push_back(_cSidechains);
	
	indent -= 20;
	height += 20;
	_cPeptide = new QCheckBox("allow peptide flex from 180°", this);
	_cPeptide->setChecked(true);
	_cPeptide->setGeometry(indent, height, CHECKBOX_WIDTH, 20);
	_cPeptide->hide();
	_allToggle.push_back(_cPeptide);

	toggleDisableOptions(0);
}

void StartScreen::toggleDisableOptions(int)
{
	bool state = _cRefine->isChecked();
	
	for (int i = 0; i < _allToggle.size(); i++)
	{
		_allToggle[i]->setEnabled(state);
	}
	
	if (!state)
	{
		return;
	}
	
	state = _cCbAngles->isChecked();
	
	for (int i = 0; i < _angleToggle.size(); i++)
	{
		_angleToggle[i]->setEnabled(state);
	}
	
	if (state)	
	{
		state = _cCgAngles->isChecked();

		for (int i = 0; i < _angle2Toggle.size(); i++)
		{
			_angle2Toggle[i]->setEnabled(state);
		}
	}

}

StartScreen::StartScreen(QWidget *parent,
                         int argc, char *argv[]) : QMainWindow(parent)
{
	this->resize(DEFAULT_WIDTH, DEFAULT_HEIGHT);

	this->setWindowTitle("Setup Vagabond");

	_argc = argc;
	_argv = argv;

	_options = OptionsPtr(new Options(_argc, (const char **)_argv));
	_options->setManual(true);
    Options::setRuntimeOptions(_options);
	_options->parse();

	_fileDialogue = NULL;

	makeButtons();
}

void StartScreen::pushRun()
{
	_options->_modelFile = _tModel->text().toStdString();
	_options->_mtzFile = _tMtz->text().toStdString();
	_options->_outputDir = _tDir->text().toStdString();
	_options->_peptideMovement = _cPeptide->isChecked();
	_options->_rSidechains = _cSidechains->isChecked();
	
	int angles = 0;
	if (_cCbAngles->isChecked())
	{
		angles = 1;
	}
	
	if (angles == 1 && _cCgAngles->isChecked())
	{
		angles = 2;
	}
	
	if (angles == 2 && _cGlyAngles->isChecked())
	{
		angles = 3;
	}
	
	_options->_bondAngles = angles;
	
	_options->_rPosition = _cPosition->isChecked();
	_options->_rInter = _cInter->isChecked();
	_options->_rIntra = _cIntra->isChecked();
	_options->_refine = _cRefine->isChecked();
	
	if (_tMaxRes->text().length())
	{
		Options::_maxRes = _tMaxRes->text().toDouble();
	}
	
	if (_tRadius->text().length())
	{
		Options::_probeRadius = _tRadius->text().toDouble();
	}
	
	if (_tRelFlex->text().length())
	{
		Options::_bStart = _tRelFlex->text().toDouble();
	}
	
	if (_tMinRes->text().length())
	{
		Options::_minRes = _tMinRes->text().toDouble();
	}
	
	VagWindow *window = new VagWindow(NULL, _argc, _argv);
	window->setStartScreen(this);
	finishUp();
	window->show();

	// memory leak of this...
}

void StartScreen::finishUp()
{
	this->hide();
//	this->deleteLater();
}

void StartScreen::getFile(std::string title, QString types,
                          QLineEdit *edit)
{
	if (_fileDialogue != NULL)
	{
		delete _fileDialogue;
	}

	_fileDialogue = new QFileDialog(this, QString::fromStdString(title),
	                                types);
	
	_fileDialogue->setNameFilter(types);

	if (types.length())
	{
		_fileDialogue->setFileMode(QFileDialog::AnyFile);
	}
	else
	{
		_fileDialogue->setFileMode(QFileDialog::DirectoryOnly);
	}

	_fileDialogue->show();

    QStringList fileNames;
    if (_fileDialogue->exec())
    {
        fileNames = _fileDialogue->selectedFiles();
    }
    
    if (fileNames.size() >= 1)
    {
		edit->setText(fileNames[0]);
    }
}

void StartScreen::chooseModel()
{
	getFile("Open model", tr("Model files (*.vbond *.pdb)"), _tModel);
}

void StartScreen::chooseMtz()
{
	getFile("Open mtz file", "Reflection files (*.mtz)", _tMtz);
}

void StartScreen::chooseDir()
{
	getFile("Open folder", "", _tDir);
}

StartScreen::~StartScreen()
{
	delete _bRun;
}

void StartScreen::toggleOptionals()
{
	_showOpts = !_showOpts;

	for (int i = 0; i < _optWidgets.size(); i++)
	{
		if (_showOpts)
		{
			_optWidgets[i]->show();
		}
		else
		{
			_optWidgets[i]->hide();
		}
	}
	
	if (_showOpts)
	{
		_bOptionals->setText("Hide optional parameters");
	}
	else
	{
		_bOptionals->setText("Show optional parameters");
	}
}

