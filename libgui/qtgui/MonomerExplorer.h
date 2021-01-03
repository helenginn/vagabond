// Vagabond
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

#ifndef __Vagabond__MonomerExplorer__
#define __Vagabond__MonomerExplorer__

#include <QtWidgets/qlabel.h>
#include <QtWidgets/qspinbox.h>
#include <QtWidgets/qlistwidget.h>
#include "../../libsrc/shared_ptrs.h"
#include "AtomListItem.h"
#include "SetterEdit.h"
#include <QtWidgets/qpushbutton.h>
#include <QtWidgets/qslider.h>
#include "VagabondGLWidget.h"
#include <map>

typedef struct
{
	ParamOptionType optionType;
	int scale;
	double value;
	int isZero;
	QLabel *lOpt;
	QLabel *lVal;
	const char *unit;
} ParamOption;

typedef std::map<QSlider *, ParamOption> OptionMap;

class MonomerExplorer : public QWidget
{
	Q_OBJECT

public:
	MonomerExplorer(QWidget *parent = NULL, MonomerPtr monomer = MonomerPtr())
	: QWidget(parent)
	{
		initialise(monomer);
	}

	MonomerPtr getMonomer()
	{
		return _monomer;
	}
	VagabondGLWidget *getKeeper()
	{
		return _keeper;
	}

	void updateCorrelation(bool force = false);
	void setKeeper(VagabondGLWidget *keeper);
	~MonomerExplorer();
private slots:
	void clickedAtomListItem();
	void pushRefineDensity();
	void pushSidechainsToEnd();
	void pushModelPosToEnd();
	void pushRefineToEnd();
	void pushSplitBond();
	void setSliderValue();
private:
	bool checkForData();
	void initialise(MonomerPtr monomer);
	void populateList();
	void makeRefinementButtons();
	void makeSlider(ParamOptionType option, int num, QString name, int min, int max, int scale, int defVal, const char *unit);
	void applyParamOptions(SamplerPtr sampled);
	Notifiable *preparePolymer();

	MonomerPtr _monomer;
	VagabondGLWidget *_keeper;

	QListWidget *_atomList;

	QLabel *_lModel;
	QLabel *_lTorsion;
	QLabel *_lKick;
	QLabel *_lWhack;
	QLabel *_lDampen;
	QLabel *_lPhi;
	QLabel *_lPsi;
	QLabel *_lRefineOpts;

	SetterEdit *_tTorsion;
	SetterEdit *_tKick;
	SetterEdit *_tWhack;
	SetterEdit *_tDampen;
	SetterEdit *_tPhi;
	SetterEdit *_tPsi;

	QPushButton *_bRefineDensity;
	QPushButton *_bRefineToEnd;
	QPushButton *_bSidechainsToEnd;
	QPushButton *_bSqueezeToEnd;
	QPushButton *_bModelPosToEnd;
	QPushButton *_bSplitBond;
	QLabel *_lCorrel;
	
	QSpinBox *_splitNumBox;
	std::vector<QWidget *> _widgets;

	OptionMap _optionMap;
	BondPtr _bond;
};


#endif 
