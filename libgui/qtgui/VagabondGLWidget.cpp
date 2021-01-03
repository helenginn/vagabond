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
#include <QCursor>
#include <cfloat>

#define PAN_SENSITIVITY 30

VagabondGLWidget::VagabondGLWidget(QWidget *obj) : SlipGL(obj)
{
	_addingWater = false;
	timer = NULL;
	_mouseButton = Qt::NoButton;
	_lastX = 0; _lastY = 0;
	_controlPressed = false;
	_shiftPressed = false;
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

void VagabondGLWidget::mouseReleaseEvent(QMouseEvent *e)
{
	if (isRefiningManually() && e->button() == Qt::RightButton)
	{
		setMouseRefine(false);
		return;
	}
	
	if (!_moving)
	{
		// this was just a click
		double prop_x = _lastX;
		double prop_y = _lastY;
		convertCoords(&prop_x, &prop_y);
		
		findAtomAtXY(prop_x, prop_y);
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
	
	if (isRefiningManually() && e->button() == Qt::RightButton)
	{
		double prop_x = e->x();
		double prop_y = e->y();

		convertCoords(&prop_x, &prop_y);
		setModelRay(prop_x, prop_y);
		return;
	}
	
	double newX = e->x();
	double xDiff = _lastX - newX;
	double newY = e->y();
	double yDiff = _lastY - newY;
	_lastX = newX;
	_lastY = newY;

	if (_mouseButton == Qt::LeftButton)
	{
		if (_controlPressed)
		{
			panned(xDiff / 100, yDiff / 100);
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

