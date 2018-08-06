#ifndef __Vagabond__MoleListItem__
#define __Vagabond__MoleListItem__

#include <QtWidgets/qlistwidget.h>
#include "../../libsrc/shared_ptrs.h"

class MoleListItem : public QListWidgetItem
{
public:
    MoleListItem(QString string, QListWidget *list, MoleculePtr mole) : QListWidgetItem(string, list)
    {
        setMole(mole);
    }

    void setMole(MoleculePtr mole)
    {
        _mole = mole;
    }

    MoleculePtr getMole()
    {
        return _mole;
    }

private:
    MoleculePtr _mole;
};

#endif

