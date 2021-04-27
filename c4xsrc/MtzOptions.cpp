// cluster4x
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

#include "MtzOptions.h"
#include "ClusterList.h"
#include <QRadioButton>
#include <QCheckBox>
#include <QLineEdit>
#include <QPushButton>
#include <QHBoxLayout>

MtzOptions::MtzOptions(QWidget *widget) : QMainWindow(widget)
{
	QWidget *window = new QWidget();
	QVBoxLayout *vout = new QVBoxLayout();
	window->setLayout(vout);
	window->setMinimumSize(400, 0);
	setWindowTitle("MTZ Options");

	{
		QRadioButton *r;
		r = new QRadioButton("Automatically load amplitudes", NULL);

		r->setChecked(true);
		r->setObjectName("r_auto");
		vout->addWidget(r);
	}

	{
		QHBoxLayout *hout = new QHBoxLayout();
		QRadioButton *r;
		r = new QRadioButton("Take column as amplitude:", NULL);
		r->setObjectName("r_1amp");
		hout->addWidget(r);

		QLineEdit *e = new QLineEdit(NULL);
		e->setPlaceholderText("MTZ column name");
		e->setObjectName("column");
		connect(e, &QLineEdit::textEdited, r, [=]() {r->setChecked(true);});
		hout->addWidget(e);
		
		vout->addLayout(hout);
	}

	{
		QHBoxLayout *hout = new QHBoxLayout();
		QRadioButton *r;
		r = new QRadioButton("Use difference between columns:", NULL);
		r->setObjectName("r_2amp");
		hout->addWidget(r);

		{
			QLineEdit *e = new QLineEdit(NULL);
			e->setPlaceholderText("1st MTZ column name");
			e->setObjectName("column+");
			connect(e, &QLineEdit::textEdited, r, [=]() {r->setChecked(true);});
			hout->addWidget(e);
		}

		{
			QLineEdit *e = new QLineEdit(NULL);
			e->setPlaceholderText("2nd MTZ column name");
			e->setObjectName("column-");
			connect(e, &QLineEdit::textEdited, r, [=]() {r->setChecked(true);});
			hout->addWidget(e);
		}
		
		vout->addLayout(hout);
	}

	{
		QHBoxLayout *hout = new QHBoxLayout();
		QCheckBox *r;
		r = new QCheckBox("Incorporate phase from column:", NULL);
		r->setObjectName("r_phase");
		hout->addWidget(r);

		{
			QLineEdit *e = new QLineEdit(NULL);
			e->setPlaceholderText("MTZ column name");
			e->setObjectName("phase");
			connect(e, &QLineEdit::textEdited, r, [=]() {r->setChecked(true);});
			hout->addWidget(e);
		}
		
		vout->addLayout(hout);
	}
	
	{
		QPushButton *p = new QPushButton("Save and continue", NULL);
		connect(p, &QPushButton::clicked, this, &MtzOptions::save);
		vout->addWidget(p);
	}
	
	setCentralWidget(window);
	show();
}

void MtzOptions::optsToScreen()
{
	std::string column = _list->valueForKey("column");
	std::string colplus = _list->valueForKey("column+");
	std::string colminus = _list->valueForKey("column-");
	
	if (column.length())
	{
		QRadioButton *r = findChild<QRadioButton *>("r_1amp");
		r->setChecked("true");
		QLineEdit *e = findChild<QLineEdit *>("column");
		e->setText(QString::fromStdString(column));
	}
	
	if (colplus.length() || colminus.length())
	{
		QRadioButton *r = findChild<QRadioButton *>("r_2amp");
		r->setChecked("true");
		QLineEdit *e = findChild<QLineEdit *>("column+");
		e->setText(QString::fromStdString(colplus));
		e = findChild<QLineEdit *>("column-");
		e->setText(QString::fromStdString(colminus));
	}

	std::string phase = _list->valueForKey("phase");
	
	if (phase.length())
	{
		QCheckBox *r = findChild<QCheckBox *>("r_phase");
		r->setChecked("true");
		QLineEdit *e = findChild<QLineEdit *>("phase");
		e->setText(QString::fromStdString(phase));
	}
}

std::string MtzOptions::lineText(std::string lineName)
{
	QLineEdit *e = findChild<QLineEdit *>(QString::fromStdString(lineName));
	std::string str = e->text().toStdString();
	return str;
}

void MtzOptions::sortPhase()
{
	QCheckBox *r = findChild<QCheckBox *>("r_phase");
	if (r->isChecked())
	{
		std::string str = lineText("phase");
		_list->addOption("phase", str);
	}
}

void MtzOptions::save()
{
	_list->clearOptions();

	QRadioButton *r = findChild<QRadioButton *>("r_auto");
	if (!r->isChecked())
	{
		QRadioButton *r = findChild<QRadioButton *>("r_1amp");
		if (r->isChecked())
		{
			std::string str = lineText("column");
			_list->addOption("column", str);
		}

		r = findChild<QRadioButton *>("r_2amp");
		if (r->isChecked())
		{
			std::string str = lineText("column+");
			_list->addOption("column+", str);
			str = lineText("column-");
			_list->addOption("column-", str);
		}
	}

	sortPhase();
	hide();
}
