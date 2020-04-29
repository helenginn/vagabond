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
#include "Query.h"


SQLInput::SQLInput(QWidget *widget) : QMainWindow(widget)
{
	_con = NULL;
	setGeometry(100, 100, 800, 660);
	show();

	setupTable();
	
	SQLCredentials *cred = new SQLCredentials(this);
	cred->setSQLParent(this);
	cred->show();
}

SQLInput::~SQLInput()
{
	if (_con != NULL)
	{
		mysql_close(_con);
	}
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

	l = new QLabel("Refinement method", this);
	l->setGeometry(left, top, 150, 40);
	l->show();
	_bin.push_back(l);

	QComboBox *c = new QComboBox(this);
	c->setGeometry(left + 150, top, 210, 40);
	c->show();
	connect(c, QOverload<int>::of(&QComboBox::currentIndexChanged),
	        this, &SQLInput::queryAltered);
	_methods = c;
	_bin.push_back(c);
	
	left += 380;

	l = new QLabel("Data reduction method", this);
	l->setGeometry(left, top, 150, 40);
	l->show();
	_bin.push_back(l);

	c = new QComboBox(this);
	c->setGeometry(left + 170, top, 210, 40);
	c->show();
	_redMeth = c;
	connect(c, QOverload<int>::of(&QComboBox::currentIndexChanged),
	        this, &SQLInput::queryAltered);
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
	connect(ch, &QCheckBox::stateChanged,
	        this, &SQLInput::queryAltered);
	_successful = ch;
	_bin.push_back(ch);

	l = new QLabel("Must be successful", this);
	l->setGeometry(left + 40, top, 140, 40);
	l->show();
	_bin.push_back(l);
	
	top += 40;
	left = 10;

	QStringList list;
	list << "load?" <<  "id" << "refinement" << "reduction"
	<< "resolution (Å)" << "R work (%)" << "R free (%)" <<
	"MTZ path" << "PDB path" << "hit" << "clustered"; 
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
	connect(b, &QPushButton::clicked, this, &SQLInput::load);
	_bin.push_back(b);
	
	left = 10;
	
	l = new QLabel("Datasets found", this);
	l->setGeometry(left, top, 400, 40);
	l->show();
	_status = l;
	_bin.push_back(l);
}

void SQLInput::outputError()
{
	const char *error = mysql_error(_con);
	QMessageBox msgBox;
	std::string str = error;
	msgBox.setText(tr("Unable to connect to database with "
	                  "your given credentials. Database says: ")
	+ QString::fromStdString(error));
	msgBox.exec();

	SQLCredentials *cred = new SQLCredentials(this);
	cred->setSQLParent(this);
	cred->show();
}

void SQLInput::connectDb(QString hostname, QString database, 
                         QString username, QString password)
{
	_con = mysql_init(NULL);

	if (_con == NULL) 
	{
		outputError();
		return;
	}

	if (mysql_real_connect(_con, hostname.toStdString().c_str(), 
	                       username.toStdString().c_str(),
	                       password.toStdString().c_str(),
	                       database.toStdString().c_str(),
	                       0, NULL, 0) == NULL) 
	{
		outputError();
		return;
	}  

	{
		Query q(_con, "SELECT DISTINCT method FROM "
		        "SARS_COV_2_Analysis_v2.Refinement;");

		_methods->clear();
		_methods->addItem("all");

		for (size_t i = 0; i < q.rowCount(); i++)
		{
			QString method = q.qValue(i, 0);
			if (method == "NULL")
			{
				method = "none";
			}

			_methods->addItem(method);
		}
	}

	{
		Query q(_con, "SELECT DISTINCT "
		        "SARS_COV_2_Analysis_v2.Data_Reduction.method FROM "
		        "SARS_COV_2_Analysis_v2.Refinement JOIN "
		        "SARS_COV_2_Analysis_v2.Data_Reduction ON "
		        "SARS_COV_2_Analysis_v2.Refinement.data_reduction_id = "
		        "SARS_COV_2_Analysis_v2.Data_Reduction.data_reduction_id;");

		_redMeth->clear();
		_redMeth->addItem("all");

		for (size_t i = 0; i < q.rowCount(); i++)
		{
			QString method = q.qValue(i, 0);
			if (method == "NULL")
			{
				method = "none";
			}

			_redMeth->addItem(method);
		}
	}

}

std::string SQLInput::constructRowQuery(bool distinctCrystals)
{
	std::string query = "SELECT ";
	
	if (!distinctCrystals)
	{
		query += "refinement_id, ";
		query += "SARS_COV_2_Analysis_v2.Refinement.method as refinement_method, ";
		query += "SARS_COV_2_Analysis_v2.Data_Reduction.method";
		query += " as reduction_method, ";
		query += "resolution_cut, rwork, rfree, ";
		query += "refinement_mtz_path, final_pdb_path";
	}
	else
	{
		query += "DISTINCT crystal_id";
	}
	
	query += " FROM SARS_COV_2_Analysis_v2.Refinement";
	
	query += " JOIN SARS_COV_2_Analysis_v2.Data_Reduction ON "
	"SARS_COV_2_Analysis_v2.Refinement.data_reduction_id = "
	"SARS_COV_2_Analysis_v2.Data_Reduction.data_reduction_id";

	query += " WHERE (";
	bool changed = false;

	if (_methods->currentIndex() > 0)
	{
		changed = true;
		std::string method = _methods->currentText().toStdString();
		query += " SARS_COV_2_Analysis_v2.Refinement.method = '";
		query += method + "' AND";
	}

	if (_redMeth->currentIndex() > 0)
	{
		changed = true;
		std::string reduct = _redMeth->currentText().toStdString();
		query += " SARS_COV_2_Analysis_v2.Data_Reduction.method = '";
		query += reduct + "' AND";
	}

	if (_successful->isChecked())
	{
		changed = true;
		query += " final_pdb_path IS NOT NULL AND";
	}
	
	query.erase(query.size() - 4, 4);
	
	if (!changed)
	{
		query.erase(query.size() - 3, 3); /* fully get rid of 'where' */
	}
	else
	{
		query += ")";
	}

	query += ";";

	return query;
}

void SQLInput::queryAltered()
{
	std::string full = constructRowQuery(false);
	std::string crystals = constructRowQuery(true);

	Query q(_con, full);
	Query cryst(_con, crystals);
	
	_results->clear();
	
	QString status = QString::number(q.rowCount()) + " datasets found "
	+ " (" + QString::number(cryst.rowCount()) + " unique crystals)";
	_status->setText(status);
	
	for (size_t i = 0; i < q.rowCount(); i++)
	{
		QTreeWidgetItem *item = new QTreeWidgetItem(_results);
		item->setCheckState(0, Qt::Checked);
		item->setText(1, q.qValue(i, 0));
		item->setText(2, q.qValue(i, 1));
		item->setText(3, q.qValue(i, 2));
		item->setText(4, q.qValue(i, 3));
		item->setText(5, q.qValue(i, 4));
		item->setText(6, q.qValue(i, 5));
		item->setText(7, q.qValue(i, 6));
		item->setText(8, q.qValue(i, 7));
	}
}

void SQLInput::load()
{
	std::cout << "Loading sets..." << std::endl;

	for (int i = 0; i < _results->topLevelItemCount(); i++)
	{
		QTreeWidgetItem *item = _results->topLevelItem(i);
		bool state = (item->checkState(0) == Qt::Checked);
		
		if (!state)
		{
			continue;
		}

		DatasetPath path;
		path.refinement_id = item->text(1).toStdString();
		path.mtz_path = item->text(7).toStdString();
		path.pdb_path = item->text(8).toStdString();
		_datasets.push_back(path);
	}

	_list->load(_datasets);
	hide();
	deleteLater();
}
