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

#include "AutoCluster.h"
#include "MtzFFT.h"
#include "ClusterList.h"
#include "MtzFile.h"
#include <hcsrc/maths.h>
#include <hcsrc/FileReader.h>
#include "Group.h"
#include <iostream>
#include <QVBoxLayout>
#include <QCheckBox>
#include <QSlider>
#include <QLabel>
#include <QPushButton>

AutoCluster::AutoCluster(Group *grp, QWidget *parent) : QMainWindow(parent)
{
	_group = grp;
	_list = NULL;
	_a = 0;
	_b = 1;
	_c = 2;
	_maxMembers = 0;
	_minPoints = 4;
	_maxDistance = 0.1;
	_biggest = -1;
	_assignAll = false;
	
	QWidget *central = new QWidget(NULL);
	QVBoxLayout *vbox = new QVBoxLayout(NULL);
	central->setLayout(vbox);

	{
		QLabel *l = new QLabel("Maximum distance between "\
		                       "cluster members:", this);
		l->setObjectName("max_distance_label");
		vbox->addWidget(l);

		QSlider *slider = new QSlider(Qt::Horizontal, NULL);
		slider->setObjectName("max_distance");
		vbox->addWidget(slider);
		slider->setMinimum(0);
		slider->setMaximum(200);
		slider->setValue(100);
		connect(slider, &QSlider::sliderReleased, 
		        this, &AutoCluster::paramChange);
	}

	{
		QLabel *l = new QLabel("Points needed for core cluster membership:",
		                       this);
		l->setObjectName("min_points_label");
		vbox->addWidget(l);

		QSlider *slider = new QSlider(Qt::Horizontal, NULL);
		slider->setObjectName("min_points");
		vbox->addWidget(slider);
		slider->setMinimum(0);
		slider->setMaximum(12);
		slider->setValue(4);
		connect(slider, &QSlider::sliderReleased, 
		        this, &AutoCluster::paramChange);
	}
	
	QCheckBox *c = new QCheckBox("Assign all remaining", this);
	c->setObjectName("remaining");
	c->setChecked(false);
	vbox->addWidget(c);
	connect(c, &QCheckBox::stateChanged, 
	        this, &AutoCluster::paramChange);
	

	QPushButton *p = new QPushButton("Convert to new groups", this);
	vbox->addWidget(p);
	connect(p, &QPushButton::clicked, this, &AutoCluster::finish);

	setCentralWidget(central);
	central->setMinimumSize(400, 0);
	show();
}

void AutoCluster::getPoints()
{
	for (size_t i = 0; i < _group->mtzCount(); i++)
	{
		DataPoint point;
		memset(&point, '\0', sizeof(DataPoint));
		point.data = _group->getMtz(i);
		point.pos.x = _group->getSvdValue(i, _a);
		point.pos.y = _group->getSvdValue(i, _a);
		point.pos.z = _group->getSvdValue(i, _b);
		vec3_mult(point.pos, _group->getDiagW(i));
		point.group = -1;
		_points.push_back(point);
	}

	std::random_shuffle(_points.begin(), _points.end());
}

void AutoCluster::assignCorePoints()
{
	int count = 0;

	for (size_t i = 0; i < _points.size(); i++)
	{
		int num = 0;
		_points[i].core = false;

		vec3 origin = _points[i].pos;
		for (size_t j = 0; j < _points.size(); j++)
		{
			if (i == j)
			{
				continue;
			}

			vec3 diff = _points[j].pos - origin;
			double l = vec3_length(diff);
			if (l < _maxDistance)
			{
				_points[i].connections.push_back(j);
				num++;
			}
		}
		
		if (num >= _minPoints)
		{
			_points[i].core = true;
			count++;
		}
	}
}

void AutoCluster::assignEdgePoints()
{
	for (size_t i = 0; i < _points.size(); i++)
	{
		if (_points[i].core)
		{
			continue;
		}

		_points[i].edge = false;
		vec3 origin = _points[i].pos;

		for (size_t j = 0; j < _points[i].connections.size(); j++)
		{
			DataPoint &other = _points[_points[i].connections[j]];
			if (!other.core)
			{
				continue;
			}

			vec3 diff = other.pos - origin;
			double l = vec3_length(diff);
			if (l < _maxDistance)
			{
				_points[i].edge = true;
				break;
			}
		}
	}
}

bool AutoCluster::findMoreBoundaries()
{
	bool more = false;
	int seed = -1;

	for (size_t i = 0; i < _points.size(); i++)
	{
		if (_points[i].group < 0 && !more)
		{
			more = true;
			seed = i;
		}
		
		if (_points[i].group > _biggest)
		{
			_biggest = _points[i].group;
		}
	}
	
	if (!more)
	{ 
		return false; 
	}
	
	_biggest++;
	int next = _biggest;
	std::vector<int> idx_to_check = _points[seed].connections;
	_points[seed].group = next;
	int count = 1;
	
	while (idx_to_check.size() > 0)
	{
		std::vector<int> next_checks;
		
		for (size_t i = 0; i < idx_to_check.size(); i++)
		{
			int idx = idx_to_check[i];
			DataPoint &p = _points[idx];
			
			if (p.group >= 0 || (!p.edge && !p.core))
			{
				continue;
			}
			
			if (p.core)
			{
				next_checks.reserve(next_checks.size() + p.connections.size());
				next_checks.insert(next_checks.end(), p.connections.begin(),
				                   p.connections.end());
			}

			p.group = next;
			count++;
		}
		
		idx_to_check = next_checks;
	}
	
	return true;
}

void AutoCluster::removeTiny()
{
	_total = _biggest;

	for (size_t n = 0; n < _biggest; n++)
	{
		int num = 0;

		for (size_t i = 0; i < _points.size(); i++)
		{
			if (_points[i].group == n)
			{
				num++;
			}
		}
		
		if (num < 3)
		{
			for (size_t i = 0; i < _points.size(); i++)
			{
				if (_points[i].group == n)
				{
					_points[i].group = -1;
				}
			}
			
			_total--;
		}
		else
		{
			std::cout << "New group with " << num << " members" << std::endl;
			if (_maxMembers < num)
			{
				_maxMembers = num;
			}
		}
	
	}
}

void AutoCluster::clear()
{
	for (size_t i = 0; i < _points.size(); i++)
	{
		_points[i].core = false;
		_points[i].edge = false;
		_points[i].group = -1;
		_points[i].connections.clear();
	}
	
	_biggest = -1;
	_total = 0;
	_maxMembers = 0;
}

void AutoCluster::assignRemaining()
{
	if (!_assignAll)
	{
		return;
	}

	for (size_t i = 0; i < _points.size(); i++)
	{
		if (_points[i].group >= 0)
		{
			continue;
		}
		
		int best_group = -1;
		double min = FLT_MAX;
		vec3 origin = _points[i].pos;
		
		for (size_t j = 0; j < _points.size(); j++)
		{
			if (i == j)
			{
				continue;
			}
			
			if (_points[j].group < 0)
			{
				continue;
			}

			vec3 diff = _points[j].pos - origin;
			double l = vec3_length(diff);
			if (l < min)
			{
				min = l;
				best_group = _points[j].group;
			}
		}
		
		if (best_group >= 0)
		{
			_points[i].group = best_group;
		}
	}

}

void AutoCluster::cluster(bool autofix)
{
	if (_group->mtzCount() < 5)
	{
		std::cout << "Not enough datasets to cluster." << std::endl;
		return;
	}
	
	if (_points.size() == 0)
	{
		getPoints();
	}
	bool acceptable = false;

	while (!acceptable)
	{
		clear();
		assignCorePoints();
		assignEdgePoints();

		while (findMoreBoundaries()) { }
		_biggest++;

		removeTiny();
		acceptable = true;

		if (autofix && _total == 1 && 
		    _maxMembers >= 0.8 * (double)_points.size())
		{
			/* thresholds too relaxed to properly segregate */
			acceptable = false;
			_maxDistance *= 0.75;
			std::cout << "Too few clusters; making thresholds "\
			"more stringent and retrying" << std::endl;
		}
	}
	
	assignRemaining();

	updateParams();
}

void AutoCluster::colour()
{
	std::cout << "Total: " << _total << " groups." << std::endl;
	int num = 0;
	for (size_t i = 0; i < _points.size(); i++)
	{
		MtzFile *file = _points[i].data->getMtzFile();
		file->setColour(0, 0, 0, 1);
	}

	for (size_t n = 0; n < _biggest; n++)
	{
		double frac = (double)num / (double)_total;
		float red = frac * 1.61 * 360;
		red = fmod(red, 360);
		float green = 70; float blue = 100;
		hsv_to_rgb(red, green, blue);

		bool found = false;
		for (size_t i = 0; i < _points.size(); i++)
		{
			if (_points[i].group == n)
			{
				MtzFile *file = _points[i].data->getMtzFile();
				file->setColour(red, green, blue, 1);
				found = true;
			}
		}
		
		if (found) { num++; }
	}
	
	if (_list)
	{
		emit _list->updateSelections();
	}
}

void AutoCluster::updateParams()
{
	QSlider *max = findChild<QSlider *>("max_distance");
	max->setValue(_maxDistance * 1000);

	{
		QLabel *lMax = findChild<QLabel *>("max_distance_label");
		QString str = "Maximum distance between cluster members: ";
		str += QString::fromStdString(f_to_str(_maxDistance, 2));
		lMax->setText(str);
	}

	QSlider *min = findChild<QSlider *>("min_points");
	min->setValue(_minPoints);

	{
		QLabel *lMin = findChild<QLabel *>("min_points_label");
		QString str = "Points needed for core cluster membership: ";
		str += QString::fromStdString(i_to_str(_minPoints));
		lMin->setText(str);
	}

	QCheckBox *remaining = findChild<QCheckBox *>("remaining");
	remaining->setChecked(_assignAll);
}

void AutoCluster::paramChange()
{
	QSlider *max = findChild<QSlider *>("max_distance");
	_maxDistance = max->value();
	_maxDistance /= 1000;

	QSlider *min = findChild<QSlider *>("min_points");
	_minPoints = min->value();

	QCheckBox *remaining = findChild<QCheckBox *>("remaining");
	_assignAll = remaining->isChecked();

	std::cout << "Maximum distance: " << _maxDistance << "; "\
	"Minimum points: " << _minPoints << std::endl;
	updateParams();

	cluster(false);
	colour();

}

void AutoCluster::finish()
{
	std::vector<MtzFFTPtr> reordered;

	for (size_t n = 0; n < _biggest; n++)
	{
		std::vector<MtzFFTPtr> list;
		for (size_t i = 0; i < _points.size(); i++)
		{
			if (_points[i].group == n)
			{
				list.push_back(_points[i].data);
				reordered.push_back(_points[i].data);
			}
		}
		
		if (list.size() == 0)
		{
			continue;
		}
		
		_list->makeGroup(list);
	}

	Group *all = _list->makeGroup(reordered);
	all->setCustomName("Re-ordered on auto-clusters");
	all->updateText();
	
	hide();
	deleteLater();
}
