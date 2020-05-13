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

#ifndef __Slip__KeeperGL
#define __Slip__KeeperGL

#include <QtWidgets/qopenglwidget.h>
#include <iostream>
#include <QtGui/qopengl.h>
#include <QtGui/qopenglfunctions.h>
#include "MtzFFT.h"
#include "Group.h"

#include <libsrc/mat4x4.h>

class GLAxis;
class QLabel;
class HKLView;
class CAlphaView;
class Plot3D;
class GLObject;
class Group;

class KeeperGL : public QOpenGLWidget, QOpenGLFunctions
{
	Q_OBJECT
	
public:
	KeeperGL(QWidget *parent);
	
	void preparePanels(int n);
	void addAxes();
	void addSVDPoints(Group *ave);
	void addUCPlot(Group *ave);
	void addHKLView(VagFFTPtr fft, double scale);
	void addCAlphaView(MtzFile *file, vec3 centre);
	void addCAlphaView(Group *ave);
	void setupCamera(void);
	
	void panned(double x, double y);
	void draggedLeftMouse(double x, double y);
	void draggedRightMouse(double x, double y);
	
	mat4x4 getProjMat()
	{
		return _projMat;
	}
	
	mat4x4 getModel()
	{
		return _model;
	}
	
	void setGroup(Group *ave);
	void saveImage(std::string filename);
	
	Group *getGroup()
	{
		return _ave;
	}
	
	Plot3D *getPlot()
	{
		return _plot;
	}
	
	CAlphaView *getCAlphaView()
	{
		return _cAlphaView;
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

	virtual void mouseMoveEvent(QMouseEvent *e);
	virtual void mousePressEvent(QMouseEvent *e);
public slots:
	
protected:
	virtual void initializeGL();
	virtual void paintGL();

private:
	void initialisePrograms();
	void updateProjection();
	void updateCamera();
	void finishCAlphaView();
	
	Group *_ave;
	Qt::MouseButton _mouseButton;
	mat4x4 _model;
	mat4x4 *_store;
	std::vector<GLAxis *> _axes;
	Plot3D *_plot;
	HKLView *_hklView;
	CAlphaView *_cAlphaView;
	QLabel *_rValues;
	bool _autoCorrect;
	
	std::vector<GLObject *> _renderMe;

	bool _controlPressed;
	bool _moving;
	double _lastX;
	double _lastY;
	double _scale = 1;

	mat4x4 _rotMat;
	mat4x4 _projMat;
	float camAlpha, camBeta, camGamma;
	float zNear, zFar;
	vec3 _centre;
	vec3 _translation;
	vec3 _transOnly;
	vec3 _totalCentroid;
};


#endif
