// Cluster4x
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

#ifndef __cluster4x__CorrelLabel__
#define __cluster4x__CorrelLabel__

#include <QLabel>

class Screen;
class MatrixView;

class CorrelLabel : public QLabel
{
Q_OBJECT
public:
	CorrelLabel(QWidget *parent, MatrixView *image, Screen *screen);

protected:
	virtual void mousePressEvent(QMouseEvent *e);
	virtual void mouseMoveEvent(QMouseEvent *e);
	virtual void mouseReleaseEvent(QMouseEvent *e);
	virtual void keyPressEvent(QKeyEvent *e);
	virtual void keyReleaseEvent(QKeyEvent *e);
	int findMtzForXY(int x, int y);

	MatrixView *_image;
	Screen *_screen;
	bool _shift;
	bool _ctrl;
	int _startFile;
	int _lastFile;
};

#endif
