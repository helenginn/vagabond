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

#include "Arbitrary.h"
#include "MyDictator.h"
#include "ClusterList.h"
#include <QVBoxLayout>
#include <QLabel>
#include <QCheckBox>
#include <QPushButton>
#include <QMessageBox>

Arbitrary::Arbitrary(QWidget *parent) : QMainWindow(parent)
{
	QWidget *window = new QWidget();
	QVBoxLayout *vout = new QVBoxLayout();
	window->setLayout(vout);
	window->setMinimumSize(400, 0);

	setWindowTitle("Arbitrary data options");

	{
		QLabel *l = new QLabel("Data processing options");
		vout->addWidget(l);
	}

	{
		QCheckBox *c = new QCheckBox("Subtract mean", NULL);
		c->setChecked(true);
		c->setObjectName("mean");
		vout->addWidget(c);
	}

	{
		QCheckBox *c = new QCheckBox("Standardise deviation", NULL);
		c->setChecked(false);
		c->setObjectName("stdev");
		vout->addWidget(c);
	}
	
	{
		QPushButton *p = new QPushButton("Load data", NULL);
		connect(p, &QPushButton::clicked, this, &Arbitrary::loadData);
		vout->addWidget(p);
	}

	setCentralWidget(window);
}

void Arbitrary::loadData()
{
	{
		QCheckBox *c = findChild<QCheckBox *>("mean");
		std::string val = (c->isChecked() ? "true" : "false");
		MyDictator::setValueForKey("average", val);
	}

	{
		QCheckBox *c = findChild<QCheckBox *>("stdev");
		std::string val = (c->isChecked() ? "true" : "false");
		MyDictator::setValueForKey("stdev", val);
	}

	try
	{
		_list->loadFromVectorList(_filename);
		hide();
		deleteLater();
	}
	catch (int e)
	{
		QMessageBox msg;
		msg.setText("Could not open file.");
		msg.exec();
		return;
	}
}
