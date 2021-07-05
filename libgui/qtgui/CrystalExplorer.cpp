//  Vagabond
//  CrystalExplorer.cpp
//
//  Created by Helen Ginn on 4/8/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "CrystalExplorer.h"
#include "../../libsrc/Crystal.h"
#include "../../libsrc/Shouter.h"
#include "../../libsrc/Motion.h"
#include "MoleListItem.h"
#include "MoleculeExplorer.h"
#include <hcsrc/FileReader.h>
#include "SetterEdit.h"
#include "VagWindow.h"
#include "../../libsrc/Polymer.h"
#include "../../libsrc/WaterNetwork.h"
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
	label->setGeometry(160, height, 400, TEXT_HEIGHT);
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
	
	if (molecule->isWaterNetwork())
	{
		height += TEXT_HEIGHT;
		QPushButton *b = new QPushButton("Refine sponges", this);
		b->setGeometry(160, height, 200, TEXT_HEIGHT);
		b->show();
		connect(b, &QPushButton::clicked,
		        [=]{ refineSponges(); });
		_widgets.push_back(b);
	}
	
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
		b->setGeometry(160, height, 160, TEXT_HEIGHT);
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

		b = new QPushButton("Reset intra-motion", this);
		b->setGeometry(160, height, 160, TEXT_HEIGHT);
		b->show();
		connect(b, &QPushButton::clicked,
		        [=]{ pushResetIntra(); });
		_widgets.push_back(b);

		b = new QPushButton("Fit intramolecular motion", this);
		b->setGeometry(380, height, 200, TEXT_HEIGHT);
		b->show();
		connect(b, &QPushButton::clicked,
		        [=]{ pushFitIntra(); });
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

		height += TEXT_HEIGHT;
		b = new QPushButton("Fit rigid body", this);
		b->setGeometry(160, height, 160, TEXT_HEIGHT);
		b->show();
		connect(b, &QPushButton::clicked,
		        [=]{ pushRigidBody(); });
		_widgets.push_back(b);

		b = new QPushButton("Write B factors to PNG", this);
		b->setGeometry(380, height, 200, TEXT_HEIGHT);
		b->show();
		connect(b, &QPushButton::clicked,
		        [=]{ pushWritePNG(); });
		_widgets.push_back(b);

		height += TEXT_HEIGHT;
		b = new QPushButton("Fit backbone", this);
		b->setGeometry(160, height, 160, TEXT_HEIGHT);
		b->show();
		connect(b, &QPushButton::clicked,
		        [=]{ pushBackbone(); });
		_widgets.push_back(b);

		b = new QPushButton("Distance matrix", this);
		b->setGeometry(380, height, 200, TEXT_HEIGHT);
		b->show();
		connect(b, &QPushButton::clicked,
		        [=]{ pushDistanceMatrix(); });
		_widgets.push_back(b);

		height += TEXT_HEIGHT;

		PolymerPtr p= ToPolymerPtr(molecule);
		AnchorPtr a = p->getAnchorModel();

		if (a->motionCount() >= 1)
		{
			label = new QLabel("Rotation scale", this);
			label->setGeometry(160, height, 200, TEXT_HEIGHT);
			_widgets.push_back(label);
			label->show();

			height += TEXT_HEIGHT;
			QSlider *s = new QSlider(Qt::Horizontal, this);
			s->setGeometry(160, height, 420, TEXT_HEIGHT);
			s->setMinimum(0);
			s->setMaximum(500);
			double now = Motion::getScale(&*a->getMotion(0));
			s->setValue(now * 100);
			s->show();
			connect(s, &QSlider::valueChanged,
			        [=]{ slideScale(s); });
			_widgets.push_back(s);
		}
	}
}

void CrystalExplorer::slideScale(QSlider *s)
{
	if (!_currMole->isPolymer())
	{
		return;
	}

	PolymerPtr p= ToPolymerPtr(_currMole);
	AnchorPtr a = p->getAnchorModel();
	
	if (a->motionCount() < 1)
	{
		return;
	}
	
	double value = s->value();
	value /= 100;

	MotionPtr mot = a->getMotion(0);
    Notifiable *notify = Options::getRuntimeOptions()->getNotify();
    notify->setObject(&*mot);
    notify->setSetter(Motion::setScale, value);
    notify->setRefreshGroup(_crystal);
    notify->setInstruction(InstructionTypeSetObjectValue);
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

void CrystalExplorer::sendObject(InstructionType type)
{
	if (!_currMole->isPolymer())
	{
		return;
	}

	PolymerPtr pol = ToPolymerPtr(_currMole);
	_vagWindow->setObject(&*pol);
	_vagWindow->pushSendInstruction(type);

}

void CrystalExplorer::pushResetMotion()
{
	sendObject(InstructionTypeResetMotion);
}

void CrystalExplorer::pushResetSides()
{
	sendObject(InstructionTypeResetSides);
}

void CrystalExplorer::pushResetIntra()
{
	sendObject(InstructionTypeResetIntra);
}

void CrystalExplorer::pushFitMotion()
{
	sendObject(InstructionTypeFitTranslation);
}

void CrystalExplorer::pushFitSides()
{
	sendObject(InstructionTypeRefineDensity);
}

void CrystalExplorer::pushDistanceMatrix()
{
	sendObject(InstructionTypeDistanceMatrix);
}

void CrystalExplorer::pushFitIntra()
{
	sendObject(InstructionTypeRefineIntramolecule);
}

void CrystalExplorer::pushRigidBody()
{
	sendObject(InstructionTypeRigidBody);
}

void CrystalExplorer::pushBackbone()
{
	sendObject(InstructionTypeRefinePosToDensity);
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

void CrystalExplorer::pushWritePNG()
{
	sendObject(InstructionTypeWritePNG);
}

void CrystalExplorer::refineSponges()
{
	if (!_currMole->isWaterNetwork())
	{
		return;
	}

	WaterNetworkPtr wat = ToWaterNetworkPtr(_currMole);
	_vagWindow->setObject(&*wat);
	_vagWindow->pushSendInstruction(InstructionTypeRefineSponges);
}

CrystalExplorer::~CrystalExplorer()
{
	clearWidgets();

	delete _moleList;
	_moleList = NULL;
	
	delete _moleExplorer;
	_moleExplorer = NULL;
}
