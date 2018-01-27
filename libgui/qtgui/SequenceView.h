#ifndef __Vagabond__SequenceView__
#define __Vagabond__SequenceView__

#include <vector>
#include <QtCore/qglobal.h>
#include <QtWidgets/qapplication.h>
#include <QtWidgets/qwidget.h>
#include <QtWidgets/qpushbutton.h>
#include "../../libsrc/shared_ptrs.h"
#include <QtWidgets/qlistwidget.h>
#include "MoleculeExplorer.h"

class ResButton;

class SequenceView : public QWidget
{
    Q_OBJECT
    
public:
    SequenceView(QWidget *parent = 0, MoleculePtr molecule = MoleculePtr())
    : QWidget(parent)
    {
        setParent(parent);
        initialise(molecule);
    }

    void setExplorer(MoleculeExplorer *me)
    {
        _explorer = me;
    }

    ~SequenceView();

private slots:
	void pushResidueButton();
	
private:
    MoleculePtr _molecule;
    MoleculeExplorer *_explorer;

    void initialise(MoleculePtr molecule);
    std::vector<ResButton *> _buttons;
};

#endif
