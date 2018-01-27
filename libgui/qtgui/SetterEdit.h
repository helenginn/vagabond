#ifndef __Vagabond__SetterEdit__
#define __Vagabond__SetterEdit__

#include "../../libsrc/shared_ptrs.h"
#include <QtWidgets/qlineedit.h>
#include "../../libsrc/RefinementStrategy.h"
#include "../../libsrc/Monomer.h"

class SetterEdit : public QLineEdit
{
    Q_OBJECT

public:
    SetterEdit(QWidget *parent = NULL) : QLineEdit(parent)
    {
        
    }

    void setMonomer(MonomerPtr monomer)
    {
        _monomer = monomer;
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
    MonomerPtr _monomer;
};


#endif
