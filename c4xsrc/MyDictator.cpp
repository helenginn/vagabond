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

#include "MyDictator.h"
#include "Screen.h"
#include "ClusterList.h"
#include <iostream>
#include <hcsrc/FileReader.h>

MyDictator::MyDictator() : Dictator()
{
	_screen = new Screen(NULL);
	_screen->show();
	_list = _screen->getList();
	_list->setDictator(this);
}

bool MyDictator::processRequest(std::string first, std::string last)
{
	if (first == "load-csv")
	{
		_list->getFromCSV(last);
	}
	else if (first == "max-res" || first == "--max-res")
	{
		double res = atof(last.c_str());
		_list->setResolution(res);
	}
	else if (first == "average")
	{
		to_lower(last);
		_properties["average"] = last;
	}
	else if (first == "transpose")
	{
		to_lower(last);
		_properties["transpose"] = last;
	}
	else if (first == "load-vectors")
	{
		_list->loadFromVectorList(last);
	}
	else if (first == "load-clusters")
	{
		if (!file_exists(last))
		{
			return true;
		}

		std::string contents = get_file_contents(last);
		_list->loadClusters(contents);
	}

	return true;
}

void MyDictator::help()
{
	std::cout << "Help statement" << std::endl;
}

void MyDictator::finished()
{
	std::cout << "Come to end of arguments" << std::endl;
	if (_list->fileCount() == 0)
	{
		_list->getFromUser();
	}

}
