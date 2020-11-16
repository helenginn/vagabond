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

#include <QTreeWidget>
#include <FileReader.h>
#include <iostream>
#include "ClusterList.h"
#include "Group.h"
#include "CrystFELInput.h"

void ClusterList::getFromStream()
{
	if (_geom.length() == 0)
	{
		std::cout << "Missing --geom=<file.geom>" << std::endl;
		exit(1);
	}
	
	CrystFELInput input(_stream, _geom, _spg, _res);
	input.setSkipAndMax(_skip, _max);
	if (_onlyLoad)
	{
		input.onlyLoad(_preload);
	}
	Group *grp = input.process();

	if (grp != NULL)
	{
		grp->setMaxResolution(_res);
		grp->setTopGroup();
		_widget->addTopLevelItem(grp);
		_clusters.push_back(grp);
		_widget->setCurrentItem(grp);
		grp->performAverage();
	}

	if (_preload.length())
	{
		std::string contents = get_file_contents(_preload);
		loadClusters(contents);
	}
}

