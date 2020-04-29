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

#include "Query.h"
#include <QMessageBox>

void finishWithError(MYSQL *sql)
{
	const char *error = mysql_error(sql);
	QMessageBox msgBox;
	std::string str = error;
	msgBox.setText(QString::fromStdString("Error: ") 
	               + QString::fromStdString(error));
	msgBox.exec();
	return;
}

Query::Query(MYSQL *sql, std::string query)
{
	_query = query;
	
	if (mysql_query(sql, query.c_str()))
	{
		finishWithError(sql);
	}

	MYSQL_RES *result = mysql_store_result(sql);

	if (result == NULL) 
	{
		finishWithError(sql);
	}

	int num_fields = mysql_num_fields(result);

	MYSQL_ROW row;

	while ((row = mysql_fetch_row(result))) 
	{ 
		for(int i = 0; i < num_fields; i++) 
		{ 
			std::string str = (row[i] ? row[i] : "NULL");
			_results.push_back(str);
		} 
	}

	mysql_free_result(result);
}
