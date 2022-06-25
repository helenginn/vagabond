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

#include "DimReduction.h"
#include "ClusterList.h"
#include "AveCSV.h"
#include "Group.h"
#include "MtzFile.h"
#include <iostream>
#include <algorithm>
#include <QVBoxLayout>
#include <QCheckBox>
#include <QSlider>
#include <QLabel>
#include <QPushButton>

DimReduction::DimReduction(Group *grp, QWidget *parent) : QMainWindow(parent)
{
	_group = grp;
	_count = grp->mtzCount();
	_list = NULL;
	_k = 8;
	_kll = 5;

	QWidget *central = new QWidget(NULL);
	QVBoxLayout *vbox = new QVBoxLayout(NULL);
	central->setLayout(vbox);

	{
		QLabel *l = new QLabel("Nearest neighbours calculated from k points:",
		                       this);
		l->setObjectName("k_value");
		vbox->addWidget(l);

		QSlider *slider = new QSlider(Qt::Horizontal, NULL);
		slider->setObjectName("min_points");
		vbox->addWidget(slider);
		slider->setMinimum(0);
		slider->setMaximum(12);
		slider->setValue(3);
		connect(slider, &QSlider::sliderReleased, 
		        this, &DimReduction::paramChange);
	}
	
	makeReducedAverage();
	getDistances();
}

void DimReduction::makeReducedAverage()
{
	_average = _group->getAveCSV();
}

void DimReduction::getDistances()
{
	int count = _group->mtzCount();
	setupMatrix(&_distances, count, count);
	setupMatrix(&_reduced, count, count);

	for (size_t i = 0; i < count - 1; i++)
	{
		for (size_t j = i + 1; j < count; j++)
		{
			double distance = 0;
			for (size_t k = 0; k < count; k++)
			{
				double one = _group->getSvdValue(i, k);
				double two = _group->getSvdValue(j, k);
				double weight = _group->getDiagW(k);
				
				double diff = weight * (one - two);
				distance += diff * diff;
			}
			
			distance = sqrt(distance) / 10;
			
			_distances.ptrs[i][j] = distance;
			_distances.ptrs[j][i] = distance;
		}
	}
}

void DimReduction::findKLLNeighbours()
{
	getReducedDistances();
	
	for (size_t i = 0; i < _count; i++)
	{
		findKLLNeighbours(i);
	}
}

void DimReduction::findKLLNeighbours(int n)
{
	std::vector<int> indices;
	std::vector<double> average;
	average.resize(_count);

	for (size_t idx = 0; idx < _count; idx++)
	{
		double dist = _reduced.ptrs[n][idx];

		if (dist != dist)
		{
			continue;
		}

		indices.push_back(idx);

		for (size_t k = 0; k < _count; k++)
		{
			double weight = _group->getDiagW(k);
			double one = _group->getSvdValue(idx, k);
			average[k] += one * weight;
		}
	}
	
	for (size_t i = 0; i < indices.size(); i++)
	{
		std::cout << indices[i] << std::endl;
	}
	
	for (size_t i = 0; i < _count; i++)
	{
		average[i] /= (double)(indices.size());
	}

	HelenCore::SVD cc;
	setupSVD(&cc, indices.size(), indices.size());

	for (size_t i_ = 1; i_ < indices.size(); i_++)
	{
		int i = indices[i_];
		for (size_t j_ = 0; j_ < i_; j_++)
		{
			int j = indices[j_];

			double dot = 0;

			for (size_t k = 0; k < _count; k++)
			{
				double weight = _group->getDiagW(k);

				double one = _group->getSvdValue(i, k);
				double two = _group->getSvdValue(j, k);

				one = weight * one - average[k];
				two = weight * two - average[k];
				
				dot += one * two;
			}

			cc.u.ptrs[i_][j_] = dot;
			cc.u.ptrs[j_][i_] = dot;
		}
	}
	
	runSVD(&cc);

	printMatrix(&cc.u);
}

void DimReduction::getReducedDistances()
{
	int count = _group->mtzCount();

	std::vector<DistIndex> vals;
	vals.reserve(count);
	
	freeMatrix(&_reduced);
	setupMatrix(&_reduced, count, count);
	
	for (size_t i = 0; i < count * count; i++)
	{
		_reduced.vals[i] = NAN;
	}

	for (size_t i = 0; i < count; i++)
	{
		_reduced.ptrs[i][i] = 0;
	}

	for (size_t i = 0; i < count; i++)
	{
		vals.clear();

		for (size_t k = 0; k < count; k++)
		{
			double dist = _distances.ptrs[i][k];
			DistIndex di = {k, dist};
			vals.push_back(di);
		}
		
		std::sort(vals.begin(), vals.end());

		for (size_t j = 0; j < _k; j++)
		{
			size_t idx = vals[j].idx;
			_reduced.ptrs[i][idx] = vals[j].dist;
			_reduced.ptrs[idx][i] = vals[j].dist;
		}
	}
}

void DimReduction::removeAverage()
{
	size_t count = _group->mtzCount();
	double sum = 0;
	double num = 0;
	for (size_t i = 0; i < count * count; i++)
	{
		double add = _reduced.vals[i];
		if (add != add)
		{
			continue;
		}

		sum += add;
		num++;
	}
	
	sum /= num;
	
	for (size_t i = 0; i < count * count; i++)
	{
		double &add = _reduced.vals[i];
		if (add != add)
		{
			continue;
		}

		add -= sum;
	}

}

void DimReduction::fillInMissingValues()
{
	int count = _group->mtzCount();

	for (size_t k = 0; k < count; k++)
	{
		for (size_t i = 0; i < count; i++)
		{
			for (size_t j = 0; j < count; j++)
			{
				double ij = _reduced.ptrs[i][j];
				double ik = _reduced.ptrs[i][k];
				double kj = _reduced.ptrs[k][j];
				
				double sum = ik + kj;
				
				if (ij != ij)
				{
					_reduced.ptrs[i][j] = sum;
					_reduced.ptrs[j][i] = sum;
				}
				else if (ij > sum)
				{
					_reduced.ptrs[i][j] = sum;
					_reduced.ptrs[j][i] = sum;
				}
			}
		}
	}
	
	for (size_t i = 0; i < count * count; i++)
	{
//		_reduced.vals[i] *= _reduced.vals[i];
		if (_reduced.vals[i] != _reduced.vals[i])
		{
//			std::cout << "!!!" << std::endl;
		}
	}
}

void DimReduction::sendToList()
{
	int count = _group->mtzCount();

	for (size_t i = 0; i < count; i++)
	{
		std::string filename = _group->getMtzFile(i)->metadata();
		for (size_t j = 0; j < count; j++)
		{
			std::string other = _group->getMtzFile(j)->metadata();
			_average->addValue(filename, other, _reduced.ptrs[i][j]);
		}
	}

	_group->useAverageType(AveComma);
}

void DimReduction::reduce()
{
	_average->clear();
	//getReducedDistances();
	findKLLNeighbours();
	
	fillInMissingValues();
	removeAverage();
	
	sendToList();
}

void DimReduction::paramChange()
{
	QSlider *min = findChild<QSlider *>("min_points");
	_k = min->value();
	updateParams();

	reduce();

}

void DimReduction::updateParams()
{

}
