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

void DropDown::init()
{
	this->resize(400, 300);

	_notice = new QLabel("Placeholder", this);
	_notice->setGeometry(50, 80, 300, 60);
	_notice->show();
	
	_combo = new QComboBox(this);
	_combo->setGeometry(100, 150, 200, 30);
	_combo->show();
}

DropDown::DropDown(LabelChoice &choice)
{
	init();
	_notice->setText("Could not find appropriate column label.\n"\
	                 "Please choose label for " + 
	                 QString::fromStdString(choice.wanted));
	
	_notice->setWordWrap(true);
	
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
}
