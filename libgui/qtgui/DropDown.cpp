// Vagabond
// Copyright (C) 2019 Helen Ginn
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

#include "DropDown.h"
#include "../../libsrc/Options.h"

void DropDown::init()
{
	this->resize(400, 250);

	_notice = new QLabel("Placeholder", this);
	_notice->setGeometry(50, 30, 300, 60);
	_notice->show();
	
	_combo = new QComboBox(this);
	_combo->setGeometry(100, 100, 200, 30);
	_combo->show();
	
	_done = new QPushButton("OK", this);
	_done->setGeometry(150, 150, 100, 30);
	_done->show();
	
	_vag = NULL;
	
	connect(_done, SIGNAL(clicked()), this, SLOT(pressOK()));
}

void DropDown::pressOK()
{
	std::cout << "Pressed OK" << std::endl;

	int chosen = _combo->currentIndex();
	std::string label = "";

	std::cout << "Chosen column: " << chosen << std::endl;
	
	if (chosen > 0)
	{
		label = _choice.availabels[chosen - 1];
	}
	
	std::cout << "Setting new label: " << label << std::endl;

	if (_choice.original == "FP")
	{
		Options::setLabF(label);
	}
	
	_vag->attemptLoadAndRun();
	hide();
	deleteLater();
}

DropDown::DropDown(LabelChoice &choice)
{
	init();
	_choice = choice;

	_notice->setText("Could not find appropriate column label.\n"\
	                 "Please choose label for " + 
	                 QString::fromStdString(choice.wanted));
	
	_notice->setWordWrap(true);
	
	_combo->addItem("...");
	
	for (int i = 0; i < choice.availabels.size(); i++)
	{
		std::string str = choice.availabels[i] + " (" + choice.types[i] + ")";
		_combo->addItem(QString::fromStdString(str));
	}
}

DropDown::DropDown()
{
	init();
}

DropDown::~DropDown()
{
	delete _notice;
	_notice = NULL;
	
	delete _combo;
	_combo = NULL;
	
	delete _done;
	_done = NULL;
}
