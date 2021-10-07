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

#ifndef __cluster4x__AutoCluster__
#define __cluster4x__AutoCluster__

#include "MtzFFTPtr.h"
#include <hcsrc/vec3.h>
#include <vector>
#include <QMainWindow>

class Group;
class ClusterList;

/** based on the DBSCAN algorithm */

class AutoCluster : public QMainWindow
{
Q_OBJECT
public:
	AutoCluster(Group *grp, QWidget *parent = NULL);
	void setAxes(int a, int b, int c)
	{
		_a = a;
		_b = b;
		_c = c;
	}

	void setList(ClusterList *list)
	{
		_list = list;
	}

	void cluster(bool autofix = true);
	
	void colour();
public slots:
	void paramChange();
	void finish();
private:
	typedef struct 
	{
		vec3 pos;
		MtzFFTPtr data;
		bool core;
		bool edge;
		int group;
		std::vector<int> connections;
	} DataPoint;

	void getPoints();
	void assignCorePoints();
	void assignEdgePoints();
	bool findMoreBoundaries();
	void removeTiny();
	void assignRemaining();
	void clear();
	void updateParams();

	/* biggest group ID */
	int _biggest;
	int _total;
	int _minPoints;
	int _maxMembers;
	double _maxDistance;

	bool _assignAll;

	int _a, _b, _c;
	Group *_group;
	ClusterList *_list;

	std::vector<DataPoint> _points;

};

#endif
