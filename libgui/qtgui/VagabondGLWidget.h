//
//  VagabondGLWidget.h
//  Vagabond
//
//  Created by Helen Ginn on 20/01/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#ifndef __Vagabond__VagabondGLWidget_h__
#define __Vagabond__VagabondGLWidget_h__

#include <h3dsrc/SlipGL.h>
#include <QtCore/qtimer.h>
#include <QtWidgets/qopenglwidget.h>
#include <QMouseEvent>
#include "VagWindow.h"
#include "../Bonds2GL.h"
#include "../Multi2GL.h"
#include "../BlobMesh.h"

class Mesh;
class Blob;
class QThread;

class VagabondGLWidget : public SlipGL
{
Q_OBJECT
public:
	VagabondGLWidget(QWidget *obj = NULL);
	~VagabondGLWidget();

	void toggleBondView();
	void toggleVisibleDensity();
	Density2GLPtr activeDensity();
	void manualRefine();
	AtomPtr findAtomAtXY(double x, double y);

	void startTimer()
	{
		if (!_timer) return;
		_timer->start();
	}

	void stopTimer()
	{
		if (!_timer) return;
		_timer->start(100);
	}

	void renderDensity(CrystalPtr crystal);
	void setDisableDensityUpdate();
	void makeMeshRightClickMenu(QPoint p);

	void setAddingWater();
	
	void setVagWindow(VagWindow *vag)
	{
		_vag = vag;
	}

	Vagabond2GLPtr getMulti2GL()
	{
		return _multi2GL;
	}
	
	Density2GLPtr getDiffDens2GL()
	{
		return _diffDens2GL;
	}
	
	Density2GLPtr getDensity2GL()
	{
		return _density2GL;
	}
	
	Density2GLPtr getWire2GL()
	{
		return _wire2GL;
	}
	
	Density2GLPtr getOrig2GL()
	{
		return _orig2GL;
	}
	
	Density2GLPtr getOrigDiff2GL()
	{
		return _origDiff2GL;
	}
	
	int getDensityState()
	{
		return _densityState;
	}
	
	Selected2GLPtr getSelectedGL()
	{
		return _selected2GL;
	}
	
	bool isRefiningManually();
	void cancelRefine();
	void resetSelection();
	void setMouseRefine(bool val);
	
	void addMesh(BlobMesh *m);
	
	void toggleKicks();
	void novalentSelected();
	void deleteSelected();
	void splitSelected();
	void focusOnSelected();
	void selectResidue(std::string chain, int number);
	void advanceMonomer(int dir);
	void setAdding(bool val);
	void setModelRay(double x, double y);
signals:
	void refine();

public slots:
	void removeBlob();
	void refineMesh();
	void smoothMesh();
	void triangulateMesh();
	void addMeshToCrystal();
	void loadBlob(Blob *b);
	void meshDensityScale(double mult);

	void meshDone();
protected:
	virtual void keyPressEvent(QKeyEvent *event);
	virtual void keyReleaseEvent(QKeyEvent *event);
	virtual void mousePressEvent(QMouseEvent *e);
	virtual void mouseReleaseEvent(QMouseEvent *e);
	virtual void mouseMoveEvent(QMouseEvent *e);

	void convertCoords(double *x, double *y);

private:
	bool prepareWorkForObject(QObject *object);

	QThread *_worker;
	Qt::MouseButton _mouseButton;
	bool _controlPressed;
	bool _shiftPressed;
	double _lastX; double _lastY;
	bool _moving;
	bool _addingWater;

	Density2GLPtr _density2GL;
	Density2GLPtr _wire2GL;
	Density2GLPtr _orig2GL;
	Density2GLPtr _origDiff2GL;
	Density2GLPtr _diffDens2GL;
	std::mutex _setup;

	int _densityState;
	
	Vagabond2GLPtr _allBond2GL;
	Vagabond2GLPtr _aveBond2GL;
	Vagabond2GLPtr _atoms2GL;
	Vagabond2GLPtr _multi2GL;
	Selected2GLPtr _selected2GL;
	
	std::vector<BlobMesh *> _meshes;
	BlobMesh *_activeBm;
	VagWindow *_vag;
};

#endif
