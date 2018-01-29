#ifndef __Vagabond__ChainMenuAction__
#define __Vagabond__ChainMenuAction__

#include <QtWidgets/qmenubar.h>
#include "../../libsrc/Molecule.h"

class ChainMenuAction : public QAction
{
public:
    ChainMenuAction(QMenu *menu, MoleculePtr molecule) : QAction(menu)
    {
        std::string actionString = "Explore " + molecule->getChainID();
        setText(QString::fromStdString(actionString));
        menu->addAction(this);

        _molecule = molecule;
    }

    MoleculePtr getMolecule()
    {
        return _molecule;
    }

private:
    MoleculePtr _molecule;
};

#endif
