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

#ifndef __cluster4x__PlotView__
#define __cluster4x__PlotView__

#include <QWidget>

typedef enum
{
	PlotSVD,
	PlotUnitCell,
} PlotType;

class Screen;
class Group;
class KeeperGL;
class AxisScroll;
class SelectionWindow;

class PlotView : public QWidget
{
Q_OBJECT
public:
	PlotView(PlotType type, QWidget *parent);
	
	~PlotView();

	KeeperGL *keeper()
	{
		return _keeper;
	}
	
	void setScreen(Screen *scr)
	{
		_scr = scr;
	}

	void setup(Group *grp);
	virtual void resizeEvent(QResizeEvent *e);
private:
	Screen *_scr;
	PlotType _type;

	KeeperGL *_keeper;
	AxisScroll *_scroll;
	SelectionWindow *_selection;
	
	std::vector<QWidget *> _bin;

};

#endif
