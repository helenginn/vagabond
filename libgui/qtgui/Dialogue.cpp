//
//  UnitCell.cpp
//  Windexing
//
//  Created by Helen Ginn on 11/12/2016.
//  Copyright (c) 2016 Helen Ginn. All rights reserved.
//

#include "Dialogue.h"
#include <QtCore/qdebug.h>
#include <string>
#include "VagWindow.h"

#define DEFAULT_WIDTH 300
#define DEFAULT_HEIGHT 140
#define BUTTON_WIDTH 100

Dialogue::Dialogue(QWidget *parent, std::string windowText,
                   std::string labelText, std::string defaultText,
                   std::string buttonText) : QMainWindow(parent)
{
	this->resize(DEFAULT_WIDTH, DEFAULT_HEIGHT);
	_window = NULL;
	_func = NULL;

	this->setWindowTitle(windowText.c_str());

	lDialogue = new QLabel(this);
	lDialogue->setText(labelText.c_str());
	lDialogue->setGeometry(50, 20, 200, 30);
	lDialogue->show();

	tDialogue = new QLineEdit(this);
	tDialogue->setPlaceholderText(defaultText.c_str());
	tDialogue->setGeometry(50, 60, 200, 30);

	bDialogue = new QPushButton(buttonText.c_str(), this);
	bDialogue->setGeometry(75, 100, 150, 30);
	bDialogue->show();

	connect(bDialogue, SIGNAL(clicked()), this, SLOT(returnClicked()));
	connect(tDialogue, SIGNAL(returnPressed()), this, SLOT(returnClicked()));
}

void Dialogue::returnClicked()
{
	QString qText = tDialogue->text();
	std::string text = qText.toStdString();

	if (_func != NULL)
	{
		double val = atof(text.c_str());
		_window->setSetter(_func, val);
		_window->setInstruction(InstructionTypeSetObjectValue);
		_window->wakeup();
		this->hide();
		_func = NULL;
		
		return;
	}

	if (_window)
	{
		_window->receiveDialogue(_type, text);
	}
}

Dialogue::~Dialogue()
{
	tDialogue->disconnect();	
	tDialogue->deleteLater();
}
