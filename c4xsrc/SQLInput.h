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

#ifndef __cluster4x__sqlinput__
#define __cluster4x__sqlinput__

#include <string>
#include <QMainWindow>
#include <QtSql>
#include <mysql/mysql.h>

typedef struct
{
	std::string mtz_path;
	std::string pdb_path;
	long refinement_id;
} DatasetPaths;

class ClusterList;
class QComboBox;
class QCheckBox;
class QTreeWidget;

class SQLInput : public QMainWindow
{
Q_OBJECT
public:
	SQLInput(QWidget *widget = NULL);
	~SQLInput();
	
	void setList(ClusterList *list)
	{
		_list = list;
	}

	void connect(QString hostname, QString database, 
	             QString username, QString password);
private:
	void outputError();
	void setupTable();
	ClusterList *_list;

	MYSQL *_con;
	QComboBox *_methods;
	QComboBox *_cutoff;
	QCheckBox *_noHits;
	QCheckBox *_noCluster;
	QCheckBox *_successful;
	QTreeWidget *_results;
	std::vector<QWidget *> _bin;
	std::vector<DatasetPaths> _datasets;
};

#endif
