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

#include "ColumnView.h"
#include "ClusterList.h"
#include <hcsrc/FileReader.h>
#include "Group.h"
#include "AveVectors.h"
#include <QVBoxLayout>
#include <QTreeWidget>

ColumnView::ColumnView(QWidget *parent, Group *ave) : QWidget(parent)
{
	_list = NULL;
	_recalculate = false;
	AverageSet *set = ave->getWorkingSet();
	_vectors = set->vectors;
	QVBoxLayout *vout = new QVBoxLayout();
	setLayout(vout);

	QTreeWidget *view = new QTreeWidget(NULL);
	vout->addWidget(view);
	
	QStringList labels;
	for (size_t i = 0; i < _vectors->titleCount(); i++)
	{
		std::string title = _vectors->title(i);
		labels.push_back(QString::fromStdString(title));
	}
	view->setHeaderLabels(labels);

	QTreeWidgetItem *item = new QTreeWidgetItem(view);
	for (size_t i = 0; i < _vectors->titleCount(); i++)
	{
		bool ch = _vectors->enabled(i);
		item->setCheckState(i, ch ? Qt::Checked : Qt::Unchecked);
	}
	view->addTopLevelItem(item);
	
	for (size_t i = 0; i < _vectors->vectorCount(); i++)
	{
		std::vector<double> v = _vectors->vector(i);
		QTreeWidgetItem *item = new QTreeWidgetItem(view);
		
		for (size_t j = 0; j < v.size(); j++)
		{
			std::string f = f_to_str(v[j], 3);
			item->setText(j, QString::fromStdString(f));
		}
		
		view->addTopLevelItem(item);
	}
	
	connect(view, &QTreeWidget::itemChanged, 
	        this, &ColumnView::itemChanged);

	show();
}

void ColumnView::itemChanged(QTreeWidgetItem *item, int column)
{
	bool checked = item->checkState(column);
	_vectors->setEnabled(column, checked);
	_recalculate = true;
}

void ColumnView::hideEvent(QHideEvent *e)
{
	if (_recalculate && _list != NULL)
	{
		_list->cluster(_list->getLastAverage());
		_recalculate = false;
	}
}
