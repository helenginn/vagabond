// Slip n Slide
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

#ifndef __cluster4x__KeeperGL__
#define __cluster4x__KeeperGL__

#include <QtWidgets/qopenglwidget.h>
#include <iostream>
#include <QtGui/qopengl.h>
#include <QtGui/qopenglfunctions.h>
#include "MtzFFT.h"
#include "Group.h"
#include <helen3d/SlipGL.h>

#include <libsrc/mat4x4.h>

class GLAxis;
class QLabel;
class HKLView;
class CAlphaView;
class ClusterPlot;
class SlipObject;
class Group;

class KeeperGL : public SlipGL
{
	Q_OBJECT
	
public:
	KeeperGL(QWidget *parent);
	
	void preparePanels(int n);
	void addAxes();
	void addPlot(Group *ave, ClusterPlot *plot);
	void addHKLView(VagFFTPtr fft, double scale);
	void addCAlphaView(MtzFile *file, vec3 centre);
	void addCAlphaView(Group *ave);
	
	void saveImage(std::string filename);
	
	ClusterPlot *getPlot()
	{
		return _plot;
	}
	
	CAlphaView *getCAlphaView()
	{
		return _cAlphaView;
	}
	
	HKLView *getHKLView()
	{
		return _hklView;
	}
	
	void setStoreMatrix(mat4x4 *store)
	{
		_model = *store;
		_store = store;
	}
	
	void setModelMatrix(mat4x4 mat)
	{
		_model = mat;
		std::cout << mat4x4_desc(mat) << std::endl;
	}

	virtual void mouseReleaseEvent(QMouseEvent *e);
	virtual void mouseMoveEvent(QMouseEvent *e);
	virtual void mousePressEvent(QMouseEvent *e);
public slots:
	
protected:

private:
	void finishCAlphaView();
	
	Qt::MouseButton _mouseButton;
	mat4x4 _model;
	mat4x4 *_store;
	std::vector<GLAxis *> _axes;
	ClusterPlot *_plot;
	HKLView *_hklView;
	CAlphaView *_cAlphaView;
	QLabel *_rValues;
	bool _autoCorrect;

	bool _controlPressed;
	bool _shiftPressed;
	double _lastX; double _lastY;
	bool _moving;
};


#endif
