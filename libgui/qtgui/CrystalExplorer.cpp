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
#include "../../libsrc/FileReader.h"
#include "SetterEdit.h"
#include "../../libsrc/Polymer.h"

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
	std::string sIntro = molecule->getClassName() + " " + molecule->getChainID();
	sIntro += " (" + i_to_str(molecule->atomCount()) + " atoms)";
	
	QString intro = QString::fromStdString(sIntro);
	int height = 0;
	
	QLabel *label = new QLabel(intro, this);
	label->setGeometry(160, height, 200, TEXT_HEIGHT);
	label->show();
	_widgets.push_back(label);
	
	std::string bfac = "Average B factor: " +
	f_to_str(molecule->getAverageBFactor(), 1);

	height += TEXT_HEIGHT;
	label = new QLabel(QString::fromStdString(bfac), this);
	label->setGeometry(160, height, 200, TEXT_HEIGHT);
	label->show();
	_widgets.push_back(label);
	
	height += TEXT_HEIGHT;
	label = new QLabel("Multiply backbone kick:", this);
	label->setGeometry(160, height, 200, TEXT_HEIGHT);
	label->show();
	_widgets.push_back(label);
	
	if (molecule->isPolymer())
	{
		SetterEdit *edit = new SetterEdit(this);
		edit->setGeometry(340, height, 100, TEXT_HEIGHT);
		edit->setText("1");

		void *parser = &*(ToParserPtr(molecule));
		edit->setSetterAndObject(parser, Polymer::vsMultiplyBackboneKick);
		edit->show();
		_widgets.push_back(edit);
	}

}

CrystalExplorer::CrystalExplorer(QWidget *parent, CrystalPtr crystal)
{
	_crystal = crystal;
	
	this->resize(600, 400);
	QString title = "Explore crystal ";
	
	title += QString::fromStdString(_crystal->getFilename());
	this->setWindowTitle(title);
	
	populateList();
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
}
