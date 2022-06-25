// vagabond
// Copyright (C) 2022 Helen Ginn
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

#ifndef __vagabond__DimReduction__
#define __vagabond__DimReduction__

#include "MtzFFTPtr.h"
#include <hcsrc/vec3.h>
#include <hcsrc/Matrix.h>
#include <vector>
#include <QMainWindow>

class Group;
class AveCSV;
class ClusterList;

class DimReduction : public QMainWindow
{
Q_OBJECT
public:
	DimReduction(Group *grp, QWidget *parent = NULL);

	void setList(ClusterList *list)
	{
		_list = list;
	}

	void reduce();
	void makeReducedAverage();
public slots:
	void paramChange();

private:
	void updateParams();
	void getDistances();
	void getReducedDistances();
	void fillInMissingValues();
	void removeAverage();
	void sendToList();
	void findKLLNeighbours();
	void findKLLNeighbours(int n);
	int _count;
	int _k;
	int _kll;

	Group *_group;
	AveCSV *_average;
	ClusterList *_list;
	
	HelenCore::Matrix _distances;
	HelenCore::Matrix _reduced;
	
	struct DistIndex
	{
		size_t idx;
		double dist;
		
		const bool operator<(const DistIndex &other) const
		{
			return (dist < other.dist);
		}
		
		const bool operator>(const DistIndex &other) const
		{
			return (dist > other.dist);
		}
	};
};

#endif
