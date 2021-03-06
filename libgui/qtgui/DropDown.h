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

#ifndef __vagabond__dropdown__
#define __vagabond__dropdown__

#include <QMainWindow>
#include <QLabel>
#include <QComboBox>
#include <QPushButton>
#include "../../libsrc/Shouter.h"
#include "VagWindow.h"

class DropDown : public QMainWindow
{
	Q_OBJECT

public:
	DropDown();
	DropDown(LabelChoice &choice);
	~DropDown();

	void setCallBack(VagWindow *vag)
	{
		_vag = vag;
	}
private slots:
	void pressOK();
	
private:
	void init();
	
	VagWindow *_vag;

	LabelChoice _choice;
	QLabel *_notice;
	QComboBox *_combo;
	QPushButton *_done;
};

#endif
