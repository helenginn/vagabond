#include <iostream>
#include "MonomerExplorer.h"
#include "../../libsrc/Monomer.h"
#include "../../libsrc/Atom.h"
#include "../../libsrc/Bond.h"
#include "../../libsrc/FileReader.h"

#define TEXT_HEIGHT 28

void MonomerExplorer::initialise(MonomerPtr monomer)
{
    _monomer = monomer;
    populateList();
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
