//  Vagabond
//
//  Created by Helen Ginn on 11/12/2016.
//  Copyright (c) 2017 Helen Ginn. All rights reserved.
//

#ifndef __Vagabond__MoleculeExplorer__
#define __Vagabond__MoleculeExplorer__

#include <stdio.h>

#include <QtCore/qglobal.h>
#include <QtWidgets/qapplication.h>
#include <QtWidgets/qwidget.h>
#include <QtWidgets/qpushbutton.h>
#include <QtWidgets/qlabel.h>
#include <QtWidgets/qtextedit.h>
#include <QtWidgets/qscrollarea.h>
#include <QtWidgets/qmainwindow.h>

#include "../../libsrc/Molecule.h"

class MonomerExplorer;
class SequenceView;
class GLKeeper;

class MoleculeExplorer : public QMainWindow
{
	Q_OBJECT

public:
	MoleculeExplorer(QWidget *parent = 0,
	                 MoleculePtr molecule = MoleculePtr());
	~MoleculeExplorer();

	void displayMonomer(MonomerPtr monomer);
	void setGLKeeper(GLKeeper *keeper);
	void updateCorrelation();

private slots:

private:
	GLKeeper *_keeper;
	MonomerExplorer *_monomerExplorer;
	MoleculePtr _molecule;
	SequenceView *_sequenceView;
	QScrollArea *_scrollArea;
};

#endif
