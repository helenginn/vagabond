//  Vagabond
//  CrystalExplorer.cpp
//
//  Created by Helen Ginn on 4/8/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "CrystalExplorer.h"
#include "../../libsrc/Crystal.h"
#include "../../libsrc/Shouter.h"
#include "MoleListItem.h"
#include "MoleculeExplorer.h"
#include "../../libsrc/FileReader.h"
#include "SetterEdit.h"
#include "VagWindow.h"
#include "../../libsrc/Polymer.h"
#include "../../libsrc/Anchor.h"

#define TEXT_HEIGHT 28

void CrystalExplorer::populateList()
{
	_moleList = new QListWidget(this);
	_moleList->setGeometry(0, 0, 150, 400);

	if (!_crystal)
	{
		shout_at_helen("Crystal object not present in "\
		               "Crystal Explorer");
		return;
	}
	
	for (int i = 0; i < _crystal->moleculeCount(); i++)
	{
		MoleculePtr mole = _crystal->molecule(i);
		QString label = QString::fromStdString(mole->getClassName());
		label += ": ";
		label += QString::fromStdString(mole->getChainID());
		new MoleListItem(label, _moleList, mole);
	}

	connect(_moleList, SIGNAL(itemSelectionChanged()), this,
	        SLOT(clickedMoleListItem()));
	_moleList->show();
}

void CrystalExplorer::clickedMoleListItem()
{
	clearWidgets();

	MoleListItem *item = static_cast<MoleListItem *>(_moleList->currentItem());
	MoleculePtr molecule = item->getMole();
	_currMole = molecule;
	std::string sIntro = molecule->getClassName() + " " + molecule->getChainID();
	sIntro += " (" + i_to_str(molecule->atomCount()) + " atoms)";
	
	QString intro = QString::fromStdString(sIntro);
	int height = 0;
	
	QLabel *label = new QLabel(intro, this);
	label->setGeometry(160, height, 200, TEXT_HEIGHT);
	label->show();
	_widgets.push_back(label);
	
	std::string bf = "Unavailable";
	
	if (!_vagWindow->isRunningSomething())
	{
		bf = f_to_str(molecule->getAverageBFactor(), 1);
	}
	
	std::string bfac = "Average B factor: " + bf;

	height += TEXT_HEIGHT;
	label = new QLabel(QString::fromStdString(bfac), this);
	label->setGeometry(160, height, 200, TEXT_HEIGHT);
	label->show();
	_widgets.push_back(label);
	
	if (molecule->isPolymer())
	{
		PolymerPtr polymer = ToPolymerPtr(molecule);
		AnchorPtr anchor = ToAnchorPtr(polymer->getAnchorModel());
		double bAnch = anchor->getBFactor();
		std::string sAnch = f_to_str(bAnch, 1);

		height += TEXT_HEIGHT;
		label = new QLabel("Anchor B factor:", this);
		label->setGeometry(160, height, 200, TEXT_HEIGHT);
		label->show();
		_widgets.push_back(label);

		SetterEdit *edit = new SetterEdit(this);
		edit->setGeometry(340, height, 100, TEXT_HEIGHT);
		edit->setText(QString::fromStdString(sAnch));
		edit->setRefreshGroup(molecule);
		edit->setSetterAndObject(&*anchor, Anchor::ssetBFactor);
		edit->show();

		_widgets.push_back(edit);
		
		QPushButton *b = new QPushButton("Sequence", this);
		b->setGeometry(460, height, 120, TEXT_HEIGHT);
		b->show();
		connect(b, &QPushButton::clicked,
		        [=]{ pushSequence(); });
		_widgets.push_back(b);
		
		height += TEXT_HEIGHT;
		b = new QPushButton("Reset motion", this);
		b->setGeometry(160, height, 120, TEXT_HEIGHT);
		b->show();
		connect(b, &QPushButton::clicked,
		        [=]{ pushResetMotion(); });
		_widgets.push_back(b);

		b = new QPushButton("Fit inter-molecular motion", this);
		b->setGeometry(380, height, 200, TEXT_HEIGHT);
		b->show();
		connect(b, &QPushButton::clicked,
		        [=]{ pushFitMotion(); });
		_widgets.push_back(b);

		height += TEXT_HEIGHT;
		b = new QPushButton("Reset sidechains", this);
		b->setGeometry(160, height, 160, TEXT_HEIGHT);
		b->show();
		connect(b, &QPushButton::clicked,
		        [=]{ pushResetSides(); });
		_widgets.push_back(b);

		b = new QPushButton("Fit sidechains", this);
		b->setGeometry(380, height, 200, TEXT_HEIGHT);
		b->show();
		connect(b, &QPushButton::clicked,
		        [=]{ pushFitSides(); });
		_widgets.push_back(b);
	}
}

CrystalExplorer::CrystalExplorer(QWidget *parent, CrystalPtr crystal)
{
	_crystal = crystal;
	_moleExplorer = NULL;
	
	this->resize(600, 400);
	QString title = "Explore crystal ";
	
	title += QString::fromStdString(_crystal->getFilename());
	this->setWindowTitle(title);
	
	populateList();
}

void CrystalExplorer::pushResetMotion()
{
	if (!_currMole->isPolymer())
	{
		return;
	}

	PolymerPtr pol = ToPolymerPtr(_currMole);
	_vagWindow->setObject(&*pol);
	_vagWindow->pushSendInstruction(InstructionTypeResetMotion);
}

void CrystalExplorer::pushResetSides()
{
	if (!_currMole->isPolymer())
	{
		return;
	}

	PolymerPtr pol = ToPolymerPtr(_currMole);
	_vagWindow->setObject(&*pol);
	_vagWindow->pushSendInstruction(InstructionTypeResetSides);
}

void CrystalExplorer::pushFitMotion()
{
	if (!_currMole->isPolymer())
	{
		return;
	}

	PolymerPtr pol = ToPolymerPtr(_currMole);
	_vagWindow->setObject(&*pol);
	_vagWindow->pushSendInstruction(InstructionTypeFitTranslation);
}

void CrystalExplorer::pushFitSides()
{
	if (!_currMole->isPolymer())
	{
		return;
	}

	PolymerPtr pol = ToPolymerPtr(_currMole);
	_vagWindow->setObject(&*pol);
	_vagWindow->pushSendInstruction(InstructionTypeRefineDensity);
}


void CrystalExplorer::pushSequence()
{
	delete _moleExplorer;
	_moleExplorer = NULL;

	_moleExplorer = new MoleculeExplorer(this, _currMole);
	_moleExplorer->setGLKeeper(_keeper);
	_moleExplorer->show();

}

void CrystalExplorer::updateCorrelation()
{
	if (_moleExplorer)
	{
		_moleExplorer->updateCorrelation();
	}
}

void CrystalExplorer::clearWidgets()
{
	for (int i = 0; i < _widgets.size(); i++)
	{
		delete _widgets[i];
		_widgets[i] = NULL;
	}
	
	_widgets.clear();
}

CrystalExplorer::~CrystalExplorer()
{
	clearWidgets();

	delete _moleList;
	_moleList = NULL;
	
	delete _moleExplorer;
	_moleExplorer = NULL;
}
