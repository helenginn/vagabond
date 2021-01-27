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

#include "PlotView.h"
#include "GLPoint.h"
#include "UCPlot.h"
#include "PlotR.h"
#include "AxisScroll.h"
#include "SelectionWindow.h"
#include "KeeperGL.h"
#include "Group.h"
#include "Screen.h"

PlotView::PlotView(PlotType type, QWidget *parent) : QWidget(parent)
{
	_type = type;
	_keeper = NULL;
	_scr = NULL;
	_scroll = NULL;
	_selection = NULL;
}

PlotView::~PlotView()
{
	for (size_t i = 0; i < _bin.size(); i++)
	{
		delete _bin[i];
	}

}

void PlotView::resizeEvent(QResizeEvent *e)
{
	int axh = 60;
	int pad = 2;

	setGeometry(0, 0, parentWidget()->width(), parentWidget()->height());

	if (_selection)
	{
		_selection->setGeometry(-pad, -pad + axh, 
		                        width() + pad, height() - axh);
	}

	if (_keeper)
	{
		_keeper->setGeometry(-pad, -pad + axh, width() + pad,
		                     height() + pad - axh);
	}

	if (_scroll)
	{
		_scroll->setGeometry(0, 0, width(), axh);
	}

}

void PlotView::setup(Group *grp)
{
	_keeper = new KeeperGL(this);
	_keeper->addAxes();
	_keeper->focusOnPosition(empty_vec3(), 20);
	
	if (_type == PlotUnitCell)
	{
		_keeper->addPlot(grp, new UCPlot());
	}
	else if (_type == PlotRFactors)
	{
		_keeper->addPlot(grp, new PlotR());
	}
	else
	{
		_keeper->addPlot(grp, new GLPoint());
	}

	Plot3D *plot = _keeper->getPlot();
	plot->repopulate();

	if (_scr)
	{
		connect(plot, &Plot3D::updateSelection,
		        _scr, &Screen::refreshSelection);
	}

	_selection = new SelectionWindow(this, _keeper);
	_selection->setPlot(plot);
	_selection->show();
	_selection->setFocusPolicy(Qt::StrongFocus);
	setFocusProxy(_selection);

	_scroll = new AxisScroll(this);
	_scroll->setPlot(plot);
	_scroll->makeLayout();
	_scroll->show();
}
