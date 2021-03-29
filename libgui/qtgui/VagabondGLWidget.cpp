//
//  VagabondGLWidget.h
//  Vagabond
//
//  Created by Helen Ginn on 20/01/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "VagabondGLWidget.h"
#include "../Density2GL.h"
#include "../Connect2GL.h"
#include "../Selected2GL.h"
#include "../../libsrc/shared_ptrs.h"
#include "../../libsrc/Options.h"
#include <QCursor>
#include "../BlobMesh.h"
#include <cfloat>
#include <QMenu>
#include <QThread>

#define PAN_SENSITIVITY 30

VagabondGLWidget::VagabondGLWidget(QWidget *obj) : SlipGL(obj)
{
	_worker = NULL;
	_addingWater = false;
	_activeBm = NULL;
	_mouseButton = Qt::NoButton;
	_lastX = 0; _lastY = 0;
	_controlPressed = false;
	_shiftPressed = false;
	zFar = 2000;
	_r = 1,
	_g = 1;
	_b = 1;

	_setup.lock();

	_centre = empty_vec3();
	_densityState = 1;

	#ifdef SETUP_BUFFERS
	setupBuffers();
	#endif // SETUP_BUFFERS

	/* Bond model render */
	_allBond2GL = Bonds2GLPtr(new Bonds2GL(false));
	
	/* Average pos render */
	_aveBond2GL = Bonds2GLPtr(new Bonds2GL(true));
	_aveBond2GL->setEnabled(false);
	
	/* Atom pos render */
	_atoms2GL = Atoms2GLPtr(new Atoms2GL());
	
	/* Atom pos render for multiple positions */
	Multi2GLPtr multi2GL = Multi2GLPtr(new Multi2GL());
	_multi2GL = multi2GL;

	/* Selected atoms render */
	_selected2GL = Selected2GLPtr(new Selected2GL());
	_selected2GL->setWidget(this);

	/* Density render */
	_density2GL = Density2GLPtr(new Density2GL());
	_density2GL->setKeeper(this);
	_density2GL->recalculate();

	/* Difference density render */
	_diffDens2GL = Density2GLPtr(new Density2GL());
	_diffDens2GL->setKeeper(this);
	_diffDens2GL->setDiffDensity(true);
	_diffDens2GL->setVisible(false);
	_diffDens2GL->recalculate();
	
	_objects.push_back(&*_allBond2GL);
	_objects.push_back(&*_aveBond2GL);
	_objects.push_back(&*_selected2GL);
	_objects.push_back(&*_atoms2GL);
	_objects.push_back(&*_multi2GL);
	_objects.push_back(&*(multi2GL->getConnected2GL()));
	_objects.push_back(&*_density2GL);
	_objects.push_back(&*_diffDens2GL);
	
	_setup.unlock();
}

void VagabondGLWidget::keyPressEvent(QKeyEvent *event)
{
	Density2GLPtr active = activeDensity();
	
	if (event->key() == Qt::Key_Alt || event->key() == Qt::Key_Control)
	{
		_controlPressed = true;
	}
	else if (event->key() == Qt::Key_Shift)
	{
		_shiftPressed = true;
		setAdding(true);
	}
	else if (event->key() == Qt::Key_D)
	{
		toggleVisibleDensity();
	}
	else if (active && (event->key() == Qt::Key_Plus ||
	                    event->key() == Qt::Key_Equal))
	{
		active->nudgeDensity(1);
	}
	else if (event->key() == Qt::Key_Minus && active)
	{
		active->nudgeDensity(-1);
	}
	else if (event->key() == Qt::Key_B)
	{
		toggleBondView();
	}
	else if (event->key() == Qt::Key_L)
	{
		_vag->toggleLog();
	}
	else if (event->key() == Qt::Key_R && !_shiftPressed)
	{
		if (!isRefiningManually())
		{
			/* start it */
			_vag->setInstructionType(InstructionTypeManualRefine);
			_vag->wakeup();
		}
		else
		{
			/* send instruction to stop asap */
			cancelRefine();
		}
	}
	else if (event->key() == Qt::Key_R && _shiftPressed)
	{
		resetSelection();
	}
	else if (event->key() == Qt::Key_Space)
	{
		focusOnSelected();
	}
	else if (event->key() == Qt::Key_S)
	{
		splitSelected();
	}
	else if (event->key() == Qt::Key_X)
	{
		deleteSelected();
	}
	else if (event->key() == Qt::Key_W && !_shiftPressed)
	{
		novalentSelected();
	}
	else if (event->key() == Qt::Key_W && _shiftPressed)
	{
		_vag->waterEverything();
	}
	else if (event->key() == Qt::Key_K)
	{
		toggleKicks();
	}
	else if (event->key() == Qt::Key_G)
	{
		_vag->gotoResidueDialogue();
	}
	else if (event->key() == Qt::Key_Comma)
	{
		advanceMonomer(-1);
	}
	else if (event->key() == Qt::Key_Period)
	{
		advanceMonomer(1);
	}
}

void VagabondGLWidget::keyReleaseEvent(QKeyEvent *event)
{
	_controlPressed = false;
	_shiftPressed = false;
	setAdding(false);
}

void VagabondGLWidget::mousePressEvent(QMouseEvent *e)
{
	_lastX = e->x();
	_lastY = e->y();
	_mouseButton = e->button();
	_moving = false;
	
	if (_addingWater && e->button() == Qt::LeftButton)
	{
		double x = e->x(); double y = e->y();
		convertCoords(&x, &y);
		setModelRay(x, y);
		bool diff = (getDensityState() == 2);
		std::cout << "Adding water on " << diff << std::endl;
		getSelectedGL()->addWater(diff);
		_addingWater = false;
		setCursor(Qt::ArrowCursor);
	}
	
	if (isRefiningManually() && e->button() == Qt::RightButton)
	{
		double x = e->x(); double y = e->y();
		convertCoords(&x, &y);
		setModelRay(x, y);
		setMouseRefine(true);
	}
}

void VagabondGLWidget::convertCoords(double *x, double *y)
{
	double w = width();
	double h = height();

	*x = 2 * *x / w - 1.0;
	*y =  - (2 * *y / h - 1.0);
}

void VagabondGLWidget::removeBlob()
{
	std::vector<BlobMesh *>::iterator it;
	it = std::find(_meshes.begin(), _meshes.end(), _activeBm);
	
	if (it != _meshes.end())
	{
		CrystalPtr cryst = Options::getActiveCrystal();
		_activeBm->deleteBlob(cryst);
		removeObject(_activeBm);
		delete *it;
		_meshes.erase(it);
		_activeBm = NULL;
	}
}

void VagabondGLWidget::makeMeshRightClickMenu(QPoint p)
{
	if (_activeBm == NULL)
	{
		return;
	}

	QMenu *menu = new QMenu;
	menu->addAction("Remove blob", this, &VagabondGLWidget::removeBlob);
	menu->addAction("Refine", this, &VagabondGLWidget::refineMesh);
	menu->addAction("Triangulate", this, &VagabondGLWidget::triangulateMesh);
	menu->addAction("Smooth", this, &VagabondGLWidget::smoothMesh);

	if (!_activeBm->hasBlob())
	{
		menu->addAction("Add to crystal", this, 
		                &VagabondGLWidget::addMeshToCrystal);
	}
	else
	{
		QMenu *scale = menu->addMenu("Scale density by...");
		scale->addAction("-50%", this, [=]() {meshDensityScale(0.5); });
		scale->addAction("-20%", this, [=]() {meshDensityScale(0.8); });
		scale->addAction("+20%", this, [=]() {meshDensityScale(1.2); });
		scale->addAction("+50%", this, [=]() {meshDensityScale(1.5); });
		
	}

	menu->exec(p);
}

void VagabondGLWidget::mouseReleaseEvent(QMouseEvent *e)
{
	if (isRefiningManually() && e->button() == Qt::RightButton)
	{
		setMouseRefine(false);
		return;
	}
	
	if (!_moving && _activeBm && e->button() == Qt::RightButton &&
	    !(_worker && _worker->isRunning()))
	{
		makeMeshRightClickMenu(e->globalPos());

		return;
	}
	
	if (!_moving)
	{
		if (_worker && _worker->isRunning())
		{
			_mouseButton = Qt::NoButton;
			return;
		}

		// this was just a click
		double prop_x = _lastX;
		double prop_y = _lastY;
		convertCoords(&prop_x, &prop_y);
		
		AtomPtr atom = findAtomAtXY(prop_x, prop_y);
		
		if (atom)
		{
			_mouseButton = Qt::NoButton;
			return;
		}
		
		double z = -FLT_MAX;
		BlobMesh *bm = NULL;
		for (size_t i = 0; i < _meshes.size(); i++)
		{
			bool found = false;
			found = _meshes[i]->intersects(prop_x, prop_y, &z);
			_meshes[i]->setSelected(false);
			
			if (found)
			{
				bm = _meshes[i];
			}

		}
		
		_activeBm = bm;
		if (bm != NULL)
		{
			bm->setSelected(true);
		}
	}
	
	_mouseButton = Qt::NoButton;
}

void VagabondGLWidget::mouseMoveEvent(QMouseEvent *e)
{
	if (_mouseButton == Qt::NoButton)
	{
		return;
	}

	_moving = true;
	
	double newX = e->x();
	double xDiff = _lastX - newX;
	double newY = e->y();
	double yDiff = _lastY - newY;
	
	if (_activeBm != NULL && !_controlPressed && 
	    !(_worker && _worker->isRunning()))
	{
		double ave = _activeBm->averageRadius();
		
		if (ave < 1)
		{
			return;
		}

		double ly = _lastY;
		double ny = newY;
		
		double diff = (ny - ly) / 10;
		
		_activeBm->setSelected(false);
		_activeBm->resize(1 + diff * 0.01);
		_activeBm->setSelectable(true);
		_activeBm->setSelected(true);

		if (_activeBm->hasBlob())
		{
			_activeBm->toBlob();
		}

		_lastX = newX;
		_lastY = newY;

		return;
	}

	_lastX = newX;
	_lastY = newY;
	
	if (isRefiningManually() && e->button() == Qt::RightButton)
	{
		double prop_x = e->x();
		double prop_y = e->y();

		convertCoords(&prop_x, &prop_y);
		setModelRay(prop_x, prop_y);
		return;
	}

	if (_mouseButton == Qt::LeftButton)
	{
		if (_controlPressed)
		{
			mat4x4 prior_model = _model;
			panned(xDiff / 100, yDiff / 100);
			vec3 empty = empty_vec3();
			double last = 1;
			
			if (_activeBm != NULL)
			{
				mat4x4 inv = mat4x4_inverse(prior_model);
				vec3 prior = mat4x4_mult_vec3(inv, empty, &last);
				updateCamera();
				inv = mat4x4_inverse(_model);
				vec3 diff = mat4x4_mult_vec3(inv, empty, &last);
				vec3_subtract_from_vec3(&diff, prior);
				_activeBm->addToVertices(diff);
				
				if (_activeBm->hasBlob())
				{
					_activeBm->toBlob();
				}

			}
		}
		else
		{
			draggedLeftMouse(xDiff * 4, yDiff * 4);
		}
	}
	else if (_mouseButton == Qt::RightButton &&
	         !isRefiningManually())
	{
		draggedRightMouse(xDiff * PAN_SENSITIVITY,
		                  yDiff * PAN_SENSITIVITY);
	}
}

void VagabondGLWidget::setDisableDensityUpdate()
{
	getDensity2GL()->setDisabled(true);
	getDiffDens2GL()->setDisabled(true);
}

void VagabondGLWidget::renderDensity(CrystalPtr crystal)
{
	getDensity2GL()->makeNewDensity(crystal);
	getDiffDens2GL()->makeNewDensity(crystal);
}

VagabondGLWidget::~VagabondGLWidget()
{

}

void VagabondGLWidget::setAddingWater()
{
	_addingWater = true;
	setCursor(Qt::CrossCursor);
}

Density2GLPtr VagabondGLWidget::activeDensity()
{
	if (_densityState == 0)
	{
		return Density2GLPtr();
	}
	else if (_densityState == 1)
	{
		return getDensity2GL();
	}
	else if (_densityState == 2)
	{
		return getDiffDens2GL();
	}
	
	return getDiffDens2GL();
}

void VagabondGLWidget::toggleVisibleDensity()
{
	if (_densityState == 0)
	{
		_densityState++;
		getDensity2GL()->setVisible(true);
		getDiffDens2GL()->setVisible(false);
	}
	else if (_densityState == 1)
	{
		_densityState++;
		getDensity2GL()->setVisible(true);
		getDiffDens2GL()->setVisible(true);
	}
	else
	{
		_densityState = 0;	
		getDensity2GL()->setVisible(false);
		getDiffDens2GL()->setVisible(false);
	}
}


void VagabondGLWidget::toggleBondView()
{
	bool enabled = _allBond2GL->isEnabled();
	_aveBond2GL->setEnabled(enabled);
	_allBond2GL->setEnabled(!enabled);
}

AtomPtr VagabondGLWidget::findAtomAtXY(double x, double y)
{
	double z = -FLT_MAX;
	AtomPtr chosen = AtomPtr();

	for (int i = 0; i < _objects.size(); i++)
	{
		SlipObject *obj = _objects[i];

		if (dynamic_cast<Vagabond2GL *>(obj) == NULL)
		{
			continue;
		}
		
		Vagabond2GL *ptr = static_cast<Vagabond2GL *>(obj);
		AtomPtr atom = ptr->findAtomAtXY(x, y, &z);
		
		if (atom)
		{
			chosen = atom;
		}
	}
	
	_selected2GL->setPicked(chosen);

	return chosen;
}

void VagabondGLWidget::manualRefine()
{
	_selected2GL->manualRefine();
}

void VagabondGLWidget::cancelRefine()
{
	_selected2GL->cancelRefine();
}

bool VagabondGLWidget::isRefiningManually()
{
	return _selected2GL->isRefining();
}

void VagabondGLWidget::setModelRay(double x, double y)
{
	/* assume a z position of -1 */
	/*
	float aspect = height / width;
	y *= aspect;
	*/
	vec3 ray = make_vec3(-x, y, -1);
	_selected2GL->setMouseRay(ray);
}

void VagabondGLWidget::setMouseRefine(bool val)
{
	_selected2GL->setMouseRefinement(val);
}

void VagabondGLWidget::focusOnSelected()
{
	_selected2GL->focusOnGroup();
}

void VagabondGLWidget::splitSelected()
{
	_selected2GL->splitSelected();
}

void VagabondGLWidget::deleteSelected()
{
	_selected2GL->deleteSelected();
}

void VagabondGLWidget::toggleKicks()
{
	_selected2GL->toggleKicks();
}

void VagabondGLWidget::advanceMonomer(int dir)
{
	_selected2GL->advanceMonomer(dir);
}

void VagabondGLWidget::setAdding(bool val)
{
	_selected2GL->setAdding(val);
}

void VagabondGLWidget::selectResidue(std::string chain, int number)
{
	_selected2GL->selectResidue(chain, number);	
}

void VagabondGLWidget::novalentSelected()
{
	_selected2GL->novalentSelected(this);
}

void VagabondGLWidget::resetSelection()
{
	_selected2GL->resetSelection();
}


void VagabondGLWidget::addMesh(BlobMesh *m)
{
	_meshes.push_back(m);
	addObject(m);
}

bool VagabondGLWidget::prepareWorkForObject(QObject *object)
{
	if (object == NULL)
	{
		return false;
	}

	if (_worker && _worker->isRunning())
	{
		std::cout << "Waiting for worker to finish old job" << std::endl;
		_worker->wait();
	}
	
	if (!_worker)
	{
		_worker = new QThread();
	}

	object->moveToThread(_worker);

	return true;
}

void VagabondGLWidget::meshDone()
{
	disconnect(this, SIGNAL(refine()), nullptr, nullptr);
	disconnect(_activeBm, SIGNAL(resultReady()), nullptr, nullptr);

	_worker->quit();
	_density2GL->allowRecalculation(true);
	_activeBm->setSelectable(true);
	if (_activeBm->hasBlob())
	{
		_activeBm->toBlob();
	}
}

void VagabondGLWidget::smoothMesh()
{
	if (_activeBm == NULL)
	{
		return;
	}

	_activeBm->setSelected(false);
	prepareWorkForObject(_activeBm);

	connect(this, SIGNAL(refine()), _activeBm, SLOT(smoothCycles()));
	connect(_activeBm, SIGNAL(resultReady()), this, SLOT(meshDone()));
	_worker->start();

	emit refine();
}

void VagabondGLWidget::triangulateMesh()
{
	if (_activeBm == NULL)
	{
		return;
	}

	_activeBm->setSelected(false);
	_activeBm->changeToTriangles();
	_activeBm->SlipObject::triangulate();
	_activeBm->changeToLines();
	_activeBm->setSelectable(true);
	_activeBm->setSelected(true);
}

void VagabondGLWidget::refineMesh()
{
	if (_activeBm == NULL)
	{
		return;
	}

	_activeBm->setSelected(false);
	_density2GL->hideUnusedVertices(_activeBm);
	_density2GL->allowRecalculation(false);
	prepareWorkForObject(_activeBm);

	connect(this, SIGNAL(refine()), _activeBm, SLOT(shrinkWrap()));
	connect(_activeBm, SIGNAL(resultReady()), this, SLOT(meshDone()));
	_worker->start();

	emit refine();
}

void VagabondGLWidget::addMeshToCrystal()
{
	if (_activeBm == NULL)
	{
		return;
	}
	
	CrystalPtr cryst = Options::getActiveCrystal();
	_activeBm->blobToCrystal(cryst);

}

void VagabondGLWidget::meshDensityScale(double mult)
{
	_activeBm->multiplyScale(mult);
}

void VagabondGLWidget::loadBlob(Blob *b)
{
	BlobMesh *bm = BlobMesh::meshFromBlob(&*_density2GL, b);
	_meshes.push_back(bm);

	_setup.lock();
	addObject(bm);
	_setup.unlock();
}
