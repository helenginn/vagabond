#include <QtWidgets/qspinbox.h>
#include <iostream>
#include "MonomerExplorer.h"
#include "../../libsrc/Monomer.h"
#include "../../libsrc/Polymer.h"
#include "../../libsrc/Atom.h"
#include "../../libsrc/Bond.h"
#include "../../libsrc/Whack.h"
#include "../../libsrc/FileReader.h"
#include "../../libsrc/Notifiable.h"
#include "../../libsrc/Options.h"

#define TEXT_HEIGHT 28

void MonomerExplorer::initialise(MonomerPtr monomer)
{
	_monomer = monomer;
	_bRefineDensity = NULL;
	_bSidechainsToEnd = NULL;
	_bRefineToEnd = NULL;
	_bModelPosToEnd = NULL;
	_bSplitBond = NULL;
	_lCorrel = NULL;
	_lModel = NULL;
	_lTorsion = NULL;
	_tTorsion = NULL;
	_lKick = NULL;
	_tKick = NULL;
	_lWhack = NULL;
	_tWhack = NULL;
	_lPhi = NULL;
	_tPhi = NULL;
	_lPsi = NULL;
	_tPsi = NULL;
	_lRefineOpts = NULL;
	_atomList = NULL;
	populateList();
	makeRefinementButtons();
}

void MonomerExplorer::setSliderValue()
{
	QObject *obj = QObject::sender();
	QSlider *slider = static_cast<QSlider *>(obj);
	ParamOption *param = &_optionMap[slider];

	param->isZero = (slider->value() == 0);

	double value = slider->value();
	value /= (double)param->scale;
	param->value = value;
	char label[100];
	sprintf(label, "%.2f%s", value, param->unit);
	QString str = label;

	param->lVal->setText(str);
}

void MonomerExplorer::makeSlider(ParamOptionType option, int num, QString name,
                                 int min, int max, int scale, int defVal, const char *unit)
{
	ParamOption param;
	param.optionType = option;

	int height = 225 + 30 * num;

	param.lOpt = new QLabel(name, this);
	param.lOpt->setGeometry(20, height, 80, 25);
	param.lOpt->show();

	// because the QSlider needs ints, all values must be
	// divided by 100 to get degrees. Hence, default is 4º.

	QSlider *slider = new QSlider(Qt::Horizontal, this);
	slider->setMinimum(min);
	slider->setMaximum(max);
	slider->setValue(defVal);
	slider->setGeometry(70, height, 100, 25);
	slider->show();
	connect(slider, SIGNAL(valueChanged(int)), this, SLOT(setSliderValue()));

	char label[100];
	double val = (double)defVal / (double)scale;
	sprintf(label, "%.2f%s", val, unit);

	param.scale = scale;
	param.isZero = (defVal == 0);
	param.value = (double)defVal / (double)scale;
	param.unit = unit;

	param.lVal = new QLabel(label, this);
	param.lVal->setGeometry(200, height, 50, 25);
	param.lVal->show(); 

	_optionMap[slider] = param;
}

void MonomerExplorer::updateCorrelation(bool force)
{
	OptionsPtr options = Options::getRuntimeOptions();
	CrystalPtr crystal = options->getActiveCrystal();
	DiffractionPtr data = options->getActiveData();
	if (!data) return;
	
	Notifiable *notify = options->getNotify();
	bool running = notify->isRunningSomething();

	if ((crystal && !running) || (crystal && force))
	{
		crystal->clearCloseCache();

		if (!_lCorrel)
		{
			delete _lCorrel;
			
			_lCorrel = new QLabel("", this); 
			_lCorrel->setGeometry(250, 200, 150, 25);
			_lCorrel->show();
		}

		double score = -_monomer->scoreWithMap(ScoreTypeCorrel, crystal, "monomer_cc.png");
		std::string scoreString = "CC (2Fo-Fc): " + f_to_str(score, 3);

		_lCorrel->setText(QString::fromStdString(scoreString));
	}
}

void MonomerExplorer::makeRefinementButtons()
{
	delete _lRefineOpts;
	_lRefineOpts = new QLabel("Refine options:", this);
	_lRefineOpts->setGeometry(20, 200, 150, 25);
	_lRefineOpts->show();

	delete _bRefineDensity;
	_bRefineDensity = new QPushButton("Refine sidechain", this);
	_bRefineDensity->setGeometry(250, 225, 150, 25);
	_bRefineDensity->show(); 
	connect(_bRefineDensity, SIGNAL(clicked()), this, SLOT(pushRefineDensity()));

	delete _bSidechainsToEnd;
	_bSidechainsToEnd = new QPushButton("Sidechains to end", this);
	_bSidechainsToEnd->setGeometry(250, 250, 150, 25);
	_bSidechainsToEnd->show(); 
	connect(_bSidechainsToEnd, SIGNAL(clicked()), this, SLOT(pushSidechainsToEnd()));

	delete _bRefineToEnd;
	_bRefineToEnd= new QPushButton("Refine to end", this);
	_bRefineToEnd->setGeometry(250, 275, 150, 25);
	_bRefineToEnd->show(); 
	connect(_bRefineToEnd, SIGNAL(clicked()), this, SLOT(pushRefineToEnd()));

	delete _bModelPosToEnd;
	_bModelPosToEnd = new QPushButton("Model pos to end", this);
	_bModelPosToEnd->setGeometry(250, 325, 150, 25);
	_bModelPosToEnd->show(); 
	connect(_bModelPosToEnd, SIGNAL(clicked()), this, SLOT(pushModelPosToEnd()));

	updateCorrelation();

	makeSlider(ParamOptionTorsion, 0, "Torsion", 0, 200, 100, 10, "º");
	makeSlider(ParamOptionKick, 1, "Kick", 0, 100, 100, 50, "");
	makeSlider(ParamOptionMagicAngles, 2, "Phi/psi", 0, 90, 1, 20, "º");
	makeSlider(ParamOptionNumBonds, 3, "Bonds", 0, 16, 1, 3, "");

}

void MonomerExplorer::populateList()
{
	delete _atomList;
	_atomList = new QListWidget(this);
	_atomList->setGeometry(0, 0, 150, 200);

	if (!_monomer)
	{
		std::cout << "Warning: no monomer!" << std::endl;
	}

	if (!_monomer->atomCount())
	{
		return;
	}

	for (int i = 0; i < _monomer->atomCount(); i++)
	{
		AtomPtr atom = _monomer->atom(i);
		QString atomName = QString::fromStdString(atom->shortDesc());
		new AtomListItem(atomName, _atomList, atom);
	}

	connect(_atomList, SIGNAL(itemSelectionChanged()), this,
	        SLOT(clickedAtomListItem()));

	_atomList->show();

}

void MonomerExplorer::setKeeper(GLKeeper *keeper)
{
	_keeper = keeper;

	/* GLKeeper focus on the atom at hand */
	AtomPtr atom = _monomer->atom(0);
	vec3 pos = atom->getAbsolutePosition();
	_keeper->focusOnPosition(pos);
}


void makeLabelAndEdit(QWidget *me, QLabel **qlabel, SetterEdit **qtext, int row,
                      QString label, QString text, bool enabled)
{ 
	delete (*qlabel);
	(*qlabel) = new QLabel(label, me);
	(*qlabel)->setGeometry(160, TEXT_HEIGHT* row, 100, TEXT_HEIGHT);
	(*qlabel)->show();

	delete (*qtext);
	(*qtext) = new SetterEdit(me);
	(*qtext)->setGeometry(270, TEXT_HEIGHT * row, 100, TEXT_HEIGHT);

	if (!enabled)
	{
		(*qtext)->setText("N/A");
		(*qtext)->setEnabled(false);
	}
	else
	{
		(*qtext)->setText(text);
	}

	MonomerExplorer *expl = static_cast<MonomerExplorer *>(me);
	(*qtext)->setRefreshGroup(expl->getMonomer());
	(*qtext)->show();
}

void MonomerExplorer::clickedAtomListItem()
{
	AtomListItem *item = static_cast<AtomListItem *>(_atomList->currentItem());
	AtomPtr atom = item->getAtom();
	QString modelType = "Model: " + QString::fromStdString(atom->getModel()->getClassName());
	
	delete _bSplitBond;
	_bSplitBond = NULL;
	
	delete _tKick; delete _lKick; delete _lWhack; delete _tWhack;
	_tKick = NULL; _lKick = NULL; _lWhack = NULL; _tWhack = NULL;

	if (atom->getModel()->isBond())
	{
		BondPtr bond = ToBondPtr(atom->getModel());
		modelType += ", " + QString::fromStdString(bond->shortDesc());        

		int labelNum = 1;
		bool enabledBond = bond->isUsingTorsion();
		bool enabledFlex = bond->getRefineFlexibility();
		double torsion = rad2deg(Bond::getTorsion(&*bond));
		QString torsionText = QString::number(torsion);

		makeLabelAndEdit(this, &_lTorsion, &_tTorsion, labelNum,
		                 "Torsion (º):", torsionText, enabledBond);        
		labelNum++;
		_tTorsion->setSetterAndObject(&*bond, Bond::setTorsion, true);

		if (!bond->getWhack())
		{
			double kick = Bond::getKick(&*bond);
			QString kickText = QString::fromStdString(f_to_str(kick, 3));

			makeLabelAndEdit(this, &_lKick, &_tKick, labelNum, "Kick:", 
			                 kickText, enabledBond & enabledFlex);
			labelNum++;
			_tKick->setSetterAndObject(&*bond, Bond::setKick);
		}
		else if (bond->getWhack())
		{
			WhackPtr whack = bond->getWhack();
			double kick = Whack::getKick(&*whack);
			QString kickText = QString::fromStdString(f_to_str(kick, 3));

			makeLabelAndEdit(this, &_lKick, &_tKick, labelNum, "Kick:", 
			                 kickText, enabledBond);
			labelNum++;
			_tKick->setSetterAndObject(&*whack, Whack::setKick);
			_tKick->setRefreshGroup(_monomer->getPolymer());

			double value = Whack::getWhack(&*whack);
			QString whackTest = QString::fromStdString(f_to_str(value, 3));

			makeLabelAndEdit(this, &_lWhack, &_tWhack, labelNum, "Whack:", 
			                 whackTest, enabledBond);
			labelNum++;
			_tWhack->setSetterAndObject(&*whack, Whack::setWhack);
			_tWhack->setRefreshGroup(_monomer->getPolymer());
		}

		double phi = rad2deg(Bond::getMagicPhi(&*bond));
		QString phiText = QString::fromStdString(f_to_str(phi, 3));
		makeLabelAndEdit(this, &_lPhi, &_tPhi, labelNum, "Phi (º)", 
		                 phiText, enabledBond & enabledFlex);
		labelNum++;
		_tPhi->setSetterAndObject(&*bond, Bond::setMagicPhi, true);

		double psi = rad2deg(Bond::getMagicPsi(&*bond));
		QString psiText = QString::fromStdString(f_to_str(psi, 3));

		makeLabelAndEdit(this, &_lPsi, &_tPsi, labelNum, "Psi (º)", 
		                 psiText, enabledBond & enabledFlex);
		labelNum++;

		_tPsi->setSetterAndObject(&*bond, Bond::setMagicPsi, true);

		// Split bond
		QLabel *label = new QLabel("For next ", this);
		label->setGeometry(420, 0, 100, 25);
		label->show();
		_widgets.push_back(label);
		label = new QLabel(" bonds", this);
		label->setGeometry(530, 0, 100, 25);
		label->show();
		_widgets.push_back(label);
		
		_splitNumBox = new QSpinBox(this);
		_splitNumBox->setGeometry(480, 0, 50, 25);
		_splitNumBox->show();
		_splitNumBox->setRange(-1, 10);

		_bSplitBond = new QPushButton("Split bond", this);
		_bSplitBond->setGeometry(420, 28, 150, 25);
		connect(_bSplitBond, SIGNAL(clicked()), this, SLOT(pushSplitBond()));
		_bSplitBond->show();
		
		_bond = ToBondPtr(atom->getModel());
	}

	delete _lModel;
	_lModel = new QLabel(modelType, this);
	_lModel->setGeometry(160, 0, 240, TEXT_HEIGHT);
	_lModel->show();
}

void MonomerExplorer::applyParamOptions(SamplerPtr sampled)
{
	sampled->clearParams();

	for (OptionMap::iterator it = _optionMap.begin(); it != _optionMap.end(); it++)
	{
		ParamOption param = it->second;
		if (param.isZero) continue;

		sampled->addParamType(param.optionType, param.value);
		
		if (param.optionType == ParamOptionTorsion)
		{
			sampled->addParamType(ParamOptionBondAngle, param.value);
		}
	}
}

Notifiable *MonomerExplorer::preparePolymer()
{
	OptionsPtr options = Options::getRuntimeOptions();
	Notifiable *notify = options->getNotify();
	PolymerPtr polymer = _monomer->getPolymer();
	applyParamOptions(polymer);
	notify->setObject(&_monomer);
	return notify;
}

bool MonomerExplorer::checkForData()
{
	OptionsPtr options = Options::getRuntimeOptions();
	DiffractionPtr data = options->getActiveData();
	if (!data)
	{
		std::cout << "No reflection file specified; cannot refine to density!" 
		<< std::endl;
		return false;
	}

	return true;	
}

void MonomerExplorer::pushSplitBond()
{
	Notifiable *notify = preparePolymer();
	notify->setObject(&*_bond);
	int num = _splitNumBox->value();
	notify->setValue(num);
	notify->setInstruction(InstructionTypeSplitBond);
}

void MonomerExplorer::pushSidechainsToEnd()
{
	if (!checkForData())
	{
		return;	
	}

	Notifiable *notify = preparePolymer();
	notify->setInstruction(InstructionTypeSidechainsToEnd);
}

void MonomerExplorer::pushModelPosToEnd()
{
	Notifiable *notify = preparePolymer();
	notify->setInstruction(InstructionTypeModelPosToEnd);
}

void MonomerExplorer::pushRefineToEnd()
{
	if (!checkForData())
	{
		return;	
	}

	Notifiable *notify = preparePolymer();
	notify->setInstruction(InstructionTypeRefineToEnd);
}

void MonomerExplorer::pushRefineDensity()
{
	if (!checkForData())
	{
		return;	
	}

	OptionsPtr options = Options::getRuntimeOptions();
	Notifiable *notify = options->getNotify();

	CrystalPtr crystal = options->getActiveCrystal();
	applyParamOptions(_monomer->getSidechain());

	_monomer->getSidechain()->setTargetRefinement(crystal, RefinementFine);
	notify->setObject(&*_monomer->getSidechain());
	notify->setGetter(AtomGroup::refine);
	notify->setRefreshGroup(_monomer);
	notify->setInstruction(InstructionTypeGetObjectValue);
}

MonomerExplorer::~MonomerExplorer()
{
	delete _atomList;
	_atomList = NULL;

	delete _lModel;
	_lModel = NULL;

	delete _lTorsion;
	_lTorsion = NULL;

	delete _tTorsion;
	_tTorsion = NULL;
	
	delete _bSplitBond;
	_bSplitBond = NULL;
}
