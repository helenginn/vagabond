// Fuck COV
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

#ifndef __fuck__cov__SelectionWindow__
#define __fuck__cov__SelectionWindow__

#include <QGraphicsView>

class KeeperGL;
class Plot3D;

class SelectionWindow : public QGraphicsView
{
Q_OBJECT
public:
	SelectionWindow(QWidget *widget, KeeperGL *keeper);
	
	void setPlot(Plot3D *plot)
	{
		_plot = plot;
	}
	
	virtual void mouseMoveEvent(QMouseEvent *e);
	virtual void mousePressEvent(QMouseEvent *e);
	virtual void mouseReleaseEvent(QMouseEvent *e);

	virtual void keyPressEvent(QKeyEvent *event);
	virtual void keyReleaseEvent(QKeyEvent *event);
	virtual void resizeEvent(QResizeEvent *e);
private:
	KeeperGL *_keeper;
	Plot3D *_plot;

	bool _selectMode;
	bool _removeMode;
	int _startX;
	int _startY;
};


#endif
