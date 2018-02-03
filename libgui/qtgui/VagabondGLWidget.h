//
//  VagabondGLWidget.h
//  Vagabond
//
//  Created by Helen Ginn on 20/01/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include <QtCore/qtimer.h>
#include <QtWidgets/qopenglwidget.h>
#include "../GLKeeper.h"
#include <Qt3DInput/qmouseevent.h>

class VagabondGLWidget : public QOpenGLWidget
{
public:
    VagabondGLWidget(QWidget *obj = NULL);
    ~VagabondGLWidget();
    
    GLKeeper *getKeeper()
    {
        return keeper;
    }
protected:
    virtual void initializeGL();
    virtual void paintGL();
    virtual void resizeGL(int w, int h);
    
    virtual void keyPressEvent(QKeyEvent *event);
    virtual void keyReleaseEvent(QKeyEvent *event);
    virtual void mousePressEvent(QMouseEvent *e);
    virtual void mouseReleaseEvent(QMouseEvent *e);
    virtual void mouseMoveEvent(QMouseEvent *e);

private:
    GLKeeper *keeper;
    QTimer *timer;
    
    Qt::MouseButton _mouseButton;
	bool _controlPressed;
    double _lastX; double _lastY;
};
