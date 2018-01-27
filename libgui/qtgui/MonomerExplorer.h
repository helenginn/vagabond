#ifndef __Vagabond__MonomerExplorer__
#define __Vagabond__MonomerExplorer__

#include <QtWidgets/qlabel.h>
#include <QtWidgets/qlistwidget.h>
#include "../../libsrc/shared_ptrs.h"
#include "AtomListItem.h"
#include "SetterEdit.h"

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

    ~MonomerExplorer();
private slots:
    void clickedAtomListItem();
private:
    void initialise(MonomerPtr monomer);
    void populateList();

    MonomerPtr _monomer;

    QListWidget *_atomList;

    QLabel *_lModel;
    QLabel *_lTorsion;
    QLabel *_lKick;
    QLabel *_lDampen;
    QLabel *_lPhi;
    QLabel *_lPsi;

    SetterEdit *_tTorsion;
    SetterEdit *_tKick;
    SetterEdit *_tDampen;
    SetterEdit *_tPhi;
    SetterEdit *_tPsi;
};


#endif 
