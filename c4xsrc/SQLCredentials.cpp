// Cluster4x
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

#include "SQLCredentials.h"
#include "SQLInput.h"
#include "credentials.h"
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>

SQLCredentials::SQLCredentials(QWidget *widget) : QMainWindow(widget)
{
	setGeometry(100, 100, 400, 300);

	int top = 10;

	QLabel *l = new QLabel("SQL credentials", this);
	l->setGeometry(40, top, 380, 40);
	l->show();
	_bin.push_back(l);
	
	top += 40;

	l = new QLabel("Hostname", this);
	l->setGeometry(40, top, 140, 40);
	l->show();
	_bin.push_back(l);

	QLineEdit *line = new QLineEdit(this);
	line->setText(DEFAULT_HOSTNAME);
	line->setGeometry(200, top, 150, 40);
	line->show();
	_hostname = line;
	_bin.push_back(line);
	top += 40;

	l = new QLabel("Database name", this);
	l->setGeometry(40, top, 140, 40);
	l->show();
	_bin.push_back(l);

	line = new QLineEdit(this);
	line->setGeometry(200, top, 150, 40);
	line->setText(DEFAULT_DATABASE);
	line->show();
	_database = line;
	_bin.push_back(line);

	top += 40;

	l = new QLabel("Username", this);
	l->setGeometry(40, top, 140, 40);
	l->show();
	_bin.push_back(l);

	line = new QLineEdit(this);
	line->setGeometry(200, top, 150, 40);
	line->setText(DEFAULT_USERNAME);
	line->show();
	_username = line;
	_bin.push_back(line);

	top += 40;

	l = new QLabel("Password", this);
	l->setGeometry(40, top, 140, 40);
	l->show();
	_bin.push_back(l);

	line = new QLineEdit(this);
	line->setGeometry(200, top, 150, 40);
	line->setEchoMode(QLineEdit::Password);
	line->show();
	_password = line;
	_password->setText(DEFAULT_PASSWORD);
	_password->setFocus();
	_bin.push_back(line);

	top += 60;

	QPushButton *push = new QPushButton("Connect", this);
	push->setGeometry(250, top, 100, 40);
	push->show();
	
	connect(push, &QPushButton::clicked, this, &SQLCredentials::connectDB);

	_bin.push_back(push);
}

void SQLCredentials::connectDB()
{
	QString hostname = _hostname->text();
	QString database = _database->text();
	QString username = _username->text();
	QString password = _password->text();

	_input->connectDb(hostname, database, username, password);
	hide();
	deleteLater();
}

SQLCredentials::~SQLCredentials()
{
	for (size_t i = 0; i < _bin.size(); i++)
	{
		delete _bin[i];
	}
}
