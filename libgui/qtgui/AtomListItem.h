#ifndef __Vagabond__AtomListItem__
#define __Vagabond__AtomListItem__

#include <QtWidgets/qlistwidget.h>
#include "../../libsrc/shared_ptrs.h"

class AtomListItem : public QListWidgetItem
{
public:
    AtomListItem(QString string, QListWidget *list, AtomPtr atom) : QListWidgetItem(string, list)
    {
        setAtom(atom);
    }

    void setAtom(AtomPtr atom)
    {
        _atom = atom;
    }

    AtomPtr getAtom()
    {
        return _atom;
    }

private:
    AtomPtr _atom;
};

#endif
