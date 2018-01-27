#ifndef __Vagabond__ResButton__
#define __Vagabond__ResButton__

#include <QtWidgets/qpushbutton.h>
#include "../../libsrc/shared_ptrs.h"

class ResButton : public QPushButton
{
    Q_OBJECT
public:
    ResButton(QWidget *parent = NULL, MonomerPtr monomer = MonomerPtr())
    : QPushButton(parent)
    {
        _monomer = monomer;
        getText();
    }

	MonomerPtr getMonomer()
	{
		return _monomer;
	}

private:
    MonomerPtr _monomer;
    void getText();
};


#endif
