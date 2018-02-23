#include <iostream>
#include "MonomerExplorer.h"
#include "../../libsrc/Monomer.h"
#include "../../libsrc/Polymer.h"
#include "../../libsrc/Atom.h"
#include "../../libsrc/Bond.h"
#include "../../libsrc/FileReader.h"
#include "../../libsrc/Notifiable.h"
#include "../../libsrc/Options.h"

#define TEXT_HEIGHT 28

void MonomerExplorer::initialise(MonomerPtr monomer)
{
    _monomer = monomer;
    _bRefineDensity = NULL;
    _bRefineToEnd = NULL;
    _lCorrel = NULL;
    _lModel = NULL;
    _lTorsion = NULL;
    _tTorsion = NULL;
    _lKick = NULL;
    _tKick = NULL;
    _lDampen = NULL;
    _tDampen = NULL;
    _lPhi = NULL;
    _tPhi = NULL;
    _lPsi = NULL;
    _tPsi = NULL;
    _lRefineOpts = NULL;
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

    int height = 225 + 25 * num;

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

    delete _bRefineToEnd;
    _bRefineToEnd= new QPushButton("Refine to end", this);
    _bRefineToEnd->setGeometry(250, 250, 150, 25);
    _bRefineToEnd->show(); 
    connect(_bRefineToEnd, SIGNAL(clicked()), this, SLOT(pushRefineToEnd()));

    OptionsPtr options = Options::getRuntimeOptions();
    CrystalPtr crystal = options->getActiveCrystal();
    Notifiable *notify = options->getNotify();
    bool running = notify->isRunningSomething();

    if (crystal && !running)
    {
        delete _lCorrel;
        double score = -_monomer->scoreWithMap(ScoreTypeCorrel, crystal);
        std::string scoreString = "CC (2Fo-Fc): " + f_to_str(score, 3);

        _lCorrel = new QLabel(QString::fromStdString(scoreString), this); 
        _lCorrel->setGeometry(250, 200, 150, 25);
        _lCorrel->show();
    }

    makeSlider(ParamOptionTorsion, 0, "Torsion", 0, 200, 1000, 100, "º");
    makeSlider(ParamOptionKick, 1, "Kick", 0, 100, 100, 50, "");
    makeSlider(ParamOptionDampen, 2, "Dampen", 0, 100, 100, 25, "");
    makeSlider(ParamOptionMagicAngles, 3, "Phi/psi", 0, 90, 1, 20, "º");
    makeSlider(ParamOptionNumBonds, 4, "Bonds", 0, 8, 1, 3, "");
}

void MonomerExplorer::populateList()
{
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
        (*qtext)->setMonomer(expl->getMonomer());
        (*qtext)->show();
}

void MonomerExplorer::clickedAtomListItem()
{
    AtomListItem *item = static_cast<AtomListItem *>(_atomList->currentItem());
    AtomPtr atom = item->getAtom();
    QString modelType = "Model: " + QString::fromStdString(atom->getModel()->getClassName());

    if (atom->getModel()->isBond())
    {
        BondPtr bond = ToBondPtr(atom->getModel());
        modelType += ", " + QString::fromStdString(bond->shortDesc());        

        bool enabledBond = bond->isUsingTorsion();
        double torsion = rad2deg(Bond::getTorsion(&*bond));
        QString torsionText = QString::number(torsion);

        double dampen = Bond::getDampening(&*bond);
        QString dampenText = QString::fromStdString(f_to_str(dampen, 3));

        makeLabelAndEdit(this, &_lTorsion, &_tTorsion, 1, "Torsion (º):",
                         torsionText, enabledBond);        
        _tTorsion->setSetterAndObject(&*bond, Bond::setTorsion, true);

        double kick = Bond::getTorsionBlur(&*bond);
        QString kickText = QString::fromStdString(f_to_str(kick, 3));

        makeLabelAndEdit(this, &_lKick, &_tKick, 2, "Kick:", kickText, enabledBond);
        _tKick->setSetterAndObject(&*bond, Bond::setTorsionBlur);

        makeLabelAndEdit(this, &_lDampen, &_tDampen, 3, "Dampen:",
                         dampenText, enabledBond);
        _tDampen->setSetterAndObject(&*bond, Bond::setDampening);

        double phi = rad2deg(Bond::getMagicPhi(&*bond));
        QString phiText = QString::fromStdString(f_to_str(phi, 3));

        makeLabelAndEdit(this, &_lPhi, &_tPhi, 4, "Phi (º)", phiText, enabledBond);
        _tPhi->setSetterAndObject(&*bond, Bond::setMagicPhi, true);

        double psi = rad2deg(Bond::getMagicPsi(&*bond));
        QString psiText = QString::fromStdString(f_to_str(psi, 3));

        makeLabelAndEdit(this, &_lPsi, &_tPsi, 5, "Psi (º)", psiText, enabledBond);

        _tPsi->setSetterAndObject(&*bond, Bond::setMagicPsi, true);
    }
    
    delete _lModel;
    _lModel = new QLabel(modelType, this);
    _lModel->setGeometry(160, 0, 240, TEXT_HEIGHT);
    _lModel->show();
}

void MonomerExplorer::applyParamOptions(SamplerPtr sampled)
{
//    sampled->clearParams();

    for (OptionMap::iterator it = _optionMap.begin(); it != _optionMap.end(); it++)
    {
        ParamOption param = it->second;
        if (param.isZero) continue;
        
//        sampled->addParamType(param.optionType, param.value);
    }
}

void MonomerExplorer::pushRefineToEnd()
{
    OptionsPtr options = Options::getRuntimeOptions();
    Notifiable *notify = options->getNotify();
    PolymerPtr polymer = _monomer->getPolymer();
    applyParamOptions(polymer);
    notify->setObject(&_monomer);
    notify->setInstruction(InstructionTypeRefineToEnd);
}

void MonomerExplorer::pushRefineDensity()
{
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
}
