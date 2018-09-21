//
//  Monomer.cpp
//  vagabond
//
//  Created by Helen Ginn on 19.05.2018
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#define DEFAULT_WIDTH (600)
#define DEFAULT_HEIGHT (450)
#define TEXTBOX_X (175)
#define LABEL_WIDTH (150)
#define TEXTBOX_WIDTH (275)

#include "StartScreen.h"
#include "../../libsrc/Options.h"
#include "VagWindow.h"
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
	
	_lInputs = new QLabel("Input files and output folders:", this);
	_lInputs->setGeometry(25, height, 300, 20);
	 
	height += 40;

	_lMtz = new QLabel("Reflection file:", this);
	_lMtz->setGeometry(25, height, LABEL_WIDTH, 20);

	std::string mtzFile = _options->_mtzFile;
	_tMtz = new QLineEdit(QString::fromStdString(mtzFile), this);
	_tMtz->setGeometry(TEXTBOX_X, height, TEXTBOX_WIDTH, 25);

	_bMtz = new QPushButton("Choose MTZ", this);
	_bMtz->setGeometry(475, height, 100, 25);
	connect(_bMtz, SIGNAL(clicked()), this, SLOT(chooseMtz()));

	height += 40;

	_lModel = new QLabel("Model file:", this);
	_lModel->setGeometry(25, height, LABEL_WIDTH, 20);

	std::string modelFile = _options->_modelFile;
	_tModel = new QLineEdit(QString::fromStdString(modelFile), this);
	_tModel->setGeometry(TEXTBOX_X, height, TEXTBOX_WIDTH, 25);
	
	_bModel = new QPushButton("Choose model", this);
	_bModel->setGeometry(475, height, 100, 25);
	connect(_bModel, SIGNAL(clicked()), this, SLOT(chooseModel()));

	height += 40;

	_lDir = new QLabel("Output folder:", this);
	_lDir->setGeometry(25, height, LABEL_WIDTH, 20);

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
	
	_lTweakables = new QLabel("Optional parameters (blank is default):", this);
	_lTweakables->setGeometry(25, height, 250, 20);
	
	height += 40;

	/*
	_lKick = new QLabel("Initial kick:", this);
	_lKick->setGeometry(25, height, LABEL_WIDTH, 20);

	_tKick = new QLineEdit("0.005", this);
	_tKick->setGeometry(TEXTBOX_X, height, 120, 25);
	
	double kick = Options::getKick();
	std::string text = f_to_str(kick, 4);
	_tKick->setText(QString::fromStdString(text));

	_lKickTip = new QLabel("Tip: reduce this to 0.002 if structure is too"\
	                       " flexible, increase to 0.01 if too rigid.", this);
	_lKickTip->setWordWrap(true);
	_lKickTip->setGeometry(310, height - 2, 270, 40);
	*/

	height += 40;

	_lMinRes = new QLabel("Min resolution (Å):", this);
	_lMinRes->setGeometry(25, height, LABEL_WIDTH, 20);

	_tMinRes = new QLineEdit("", this);
	_tMinRes->setGeometry(TEXTBOX_X, height, 120, 25);
	
	if (Options::minRes() > 0)
	{
		double res = Options::minRes();
		std::string text = f_to_str(res, 3);
		_tMinRes->setText(QString::fromStdString(text));
	}

	height += 40;

	_lMaxRes = new QLabel("Max resolution (Å):", this);
	_lMaxRes->setGeometry(25, height, LABEL_WIDTH, 20);

	_tMaxRes = new QLineEdit("", this);
	_tMaxRes->setGeometry(TEXTBOX_X, height, 120, 25);

	if (Options::maxRes() > 0)
	{
		double res = Options::maxRes();
		std::string text = f_to_str(res, 3);
		_tMaxRes->setText(QString::fromStdString(text));
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
	
	if (false && _tKick->text().length())
	{
		Options::_kick = _tKick->text().toDouble();
	}
	
	if (_tMaxRes->text().length())
	{
		Options::_maxRes = _tMaxRes->text().toDouble();
	}
	
	if (_tMinRes->text().length())
	{
		Options::_minRes = _tMinRes->text().toDouble();
	}
	
	VagWindow *window = new VagWindow(NULL, _argc, _argv);
	window->show();

	this->hide();
	// memory leak of this... 
}

void StartScreen::getFile(std::string title, std::string types,
                          QLineEdit *edit)
{
	if (_fileDialogue != NULL)
	{
		delete _fileDialogue;
	}

	_fileDialogue = new QFileDialog(this, QString::fromStdString(title),
	                                QString::fromStdString(types));

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
	getFile("Open model", "Model files (*.pdb, *.vbond)", _tModel);
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
