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

#ifndef __cluster4x__query__
#define __cluster4x__query__

#include <mysql.h>
#include <string>
#include <vector>
#include <QString>

typedef std::vector<std::string> Results;

class Query
{
public:
	Query(MYSQL *con, std::string query);

	size_t rowCount()
	{
		return _rows.size();
	}

	size_t valueCount()
	{
		return _num_fields;
	}

	QString qValue(int i, int j)
	{
		return QString::fromStdString(_rows[i][j]);
	}
private:
	std::string _query;
	size_t _num_fields;

	std::vector<Results> _rows;
};

#endif
