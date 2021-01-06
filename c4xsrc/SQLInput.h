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

#include "DatasetPath.h"
#include <string>
#include <QMainWindow>
#include <QtSql>
#include <mysql/mysql.h>

class QComboBox;
class QCheckBox;
class QLabel;
class QTreeWidget;
class QSlider;
class QLineEdit;
class ClusterList;

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
	
	void setTargetID(std::string target)
	{
		_targetID = target;
	}

	void connectDb(QString hostname, QString database, 
	               QString username, QString password);
private slots:
	void queryAltered();
	void load();
	void sliderChanged();
private:
	std::string constructRowQuery(bool distinctCrystals);
	void outputError();
	void setupTable();
	void updateSlider();
	ClusterList *_list;

	MYSQL *_con;
	QComboBox *_methods;
	QComboBox *_redMeth;
	QCheckBox *_noHits;
	QCheckBox *_noCluster;
	QCheckBox *_successful;
	QSlider *_slider;
	QLineEdit *_chunker;
	QLabel *_status;
	QLabel *_chosen;
	QTreeWidget *_results;
	std::vector<QWidget *> _bin;
	std::vector<DatasetPath> _datasets;
	std::string _targetID;
	
	size_t _chunkSize;
};

#endif
