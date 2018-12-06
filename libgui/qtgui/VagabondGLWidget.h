//
//  VagabondGLWidget.h
//  Vagabond
//
//  Created by Helen Ginn on 20/01/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#ifndef __Vagabond__VagabondGLWidget_h__
#define __Vagabond__VagabondGLWidget_h__

#include <QtCore/qtimer.h>
#include <QtWidgets/qopenglwidget.h>
#include "../GLKeeper.h"
#include <Qt3DInput/qmouseevent.h>
#include "VagWindow.h"

class VagabondGLWidget : public QOpenGLWidget
{
public:
	VagabondGLWidget(QWidget *obj = NULL);
	~VagabondGLWidget();

	GLKeeper *getKeeper()
	{
		return keeper;
	}

	void startTimer()
	{
		if (!timer) return;
		timer->start();
	}

	void stopTimer()
	{
		if (!timer) return;
		timer->start(100);
	}
	virtual void resizeGL();

	void renderDensity(CrystalPtr crystal);
	void manualRefine();
	
	void setVagWindow(VagWindow *vag)
	{
		_vag = vag;
	}
protected:
	virtual void initializeGL();
	virtual void paintGL();

	virtual void keyPressEvent(QKeyEvent *event);
	virtual void keyReleaseEvent(QKeyEvent *event);
	virtual void mousePressEvent(QMouseEvent *e);
	virtual void mouseReleaseEvent(QMouseEvent *e);
	virtual void mouseMoveEvent(QMouseEvent *e);

	void convertCoords(double *x, double *y);

	private:
	GLKeeper *keeper;
	QTimer *timer;

	Qt::MouseButton _mouseButton;
	bool _controlPressed;
	double _lastX; double _lastY;
	bool _moving;
	
	VagWindow *_vag;
};

#endif
