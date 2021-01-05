#ifndef __Vagabond__SetterEdit__
#define __Vagabond__SetterEdit__

#include "../../libsrc/shared_ptrs.h"
#include <QtWidgets/qlineedit.h>
#include <hcsrc/RefinementStrategy.h>
#include "../../libsrc/Monomer.h"

class SetterEdit : public QLineEdit
{
    Q_OBJECT

public:
    SetterEdit(QWidget *parent = NULL) : QLineEdit(parent)
    {
        
    }

    void setRefreshGroup(AtomGroupPtr monomer)
    {
        _group = monomer;
    }

    void setSetterAndObject(void *object, Setter setter, bool degrees = false)
    {
        _setter = setter;
        _object = object;
        _degrees = degrees;
        connect(this, SIGNAL(editingFinished()),
                this, SLOT(setObjectValue()));
    }

private slots:
    void setObjectValue();

private:
    Setter _setter;
    void *_object;
    bool _degrees;
    AtomGroupPtr _group;
};


#endif
