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

#include <iostream>
#include <QMessageBox>
#include <QComboBox>
#include <QCheckBox>
#include <QLabel>
#include <QPushButton>
#include <QTreeWidget>

#include "ClusterList.h"
#include "SQLInput.h"
#include "SQLCredentials.h"

SQLInput::SQLInput(QWidget *widget) : QMainWindow(widget)
{
	_db = QSqlDatabase::addDatabase("QMYSQL");
	setGeometry(100, 100, 800, 660);
	show();

	setupTable();
	
	SQLCredentials *cred = new SQLCredentials(this);
	cred->setSQLParent(this);
	cred->show();
}

void SQLInput::setupTable()
{
	int left = 10;
	int top = 10;
	QLabel *l = new QLabel("Filters", this);
	l->setGeometry(left, top, 380, 40);
	l->show();
	_bin.push_back(l);

	top += 40;

	l = new QLabel("Method", this);
	l->setGeometry(left, top, 100, 40);
	l->show();
	_bin.push_back(l);

	QComboBox *c = new QComboBox(this);
	c->setGeometry(left + 100, top, 150, 40);
	c->show();
	_methods = c;
	_bin.push_back(c);
	
	left += 360;

	l = new QLabel("Res cutoff (Å)", this);
	l->setGeometry(left, top, 100, 40);
	l->show();
	_bin.push_back(l);

	c = new QComboBox(this);
	c->setGeometry(left + 100, top, 150, 40);
	c->show();
	_cutoff = c;
	_bin.push_back(c);
	
	left = 10;
	top += 40;

	QCheckBox *ch = new QCheckBox(this);
	ch->setChecked(true);
	ch->setGeometry(left, top, 40, 40);
	ch->show();
	_noHits = ch;
	_bin.push_back(ch);

	l = new QLabel("No PanDDA hits yet", this);
	l->setGeometry(left + 40, top, 160, 40);
	l->show();
	_bin.push_back(l);
	
	left += 200;

	ch = new QCheckBox(this);
	ch->setChecked(false);
	ch->setGeometry(left, top, 40, 40);
	ch->show();
	_noCluster = ch;
	_bin.push_back(ch);

	l = new QLabel("Not yet member of cluster", this);
	l->setGeometry(left + 40, top, 200, 40);
	l->show();
	_bin.push_back(l);
	
	left += 290;

	ch = new QCheckBox(this);
	ch->setChecked(true);
	ch->setGeometry(left, top, 40, 40);
	ch->show();
	_successful = ch;
	_bin.push_back(ch);

	l = new QLabel("Successful", this);
	l->setGeometry(left + 40, top, 140, 40);
	l->show();
	_bin.push_back(l);
	
	top += 40;
	left = 10;

	QStringList list;
	list << "load?" <<  "id" << "method" << "res cutoff (Å)"
	<< "resolution (Å)" << "R work (%)" << "R free (%)" << "hit" << 
	"clustered"; 
	_results = new QTreeWidget(this);
	_results->setGeometry(left, top, width() - 20, height() - top - 60);
	_results->setHeaderLabels(list);
	
	_results->resizeColumnToContents(0);
	_results->resizeColumnToContents(7);
	_results->resizeColumnToContents(8);

	_results->show();
	_bin.push_back(_results);
	
	top = height() - 50;
	left = width() - 120;
	
	QPushButton *b = new QPushButton("Load", this);
	b->setGeometry(left, top, 100, 40);
	b->show();
	_bin.push_back(b);
}

void SQLInput::connect(QString hostname, QString database, 
                       QString username, QString password)
{
	_db.setHostName(hostname);
	_db.setDatabaseName(database);
	_db.setUserName(username);
	_db.setPassword(password);
	bool ok = _db.open();

	if (!ok)
	{
		QSqlError err = _db.lastError();
		QMessageBox msgBox;
		msgBox.setText(tr("Unable to connect to database with "
		                  "your given credentials. Database says: ")
		                  + err.databaseText());
		msgBox.exec();
		
		SQLCredentials *cred = new SQLCredentials(this);
		cred->setSQLParent(this);
		cred->show();
		return;
	}
	
	QSqlQuery query = _db.exec("SELECT DISTINCT method FROM "
	                           "SARS_COV_2_Analysis_v2.Refinement;");

	_methods->clear();
	_cutoff->clear();
	_methods->addItem("all");
	_cutoff->addItem("any");

	while (query.next())
	{
		QString method = query.value(0).toString();
		if (method == "__NULL__")
		{
			method = "none";
		}

		_methods->addItem(method);
	}
}
