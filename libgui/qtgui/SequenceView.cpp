#include "SequenceView.h"
#include "ResButton.h"
#include "../../libsrc/Polymer.h"

#define BUTTON_DIM 30

void SequenceView::initialise(MoleculePtr molecule)
{
    _molecule = molecule;
    int currX = 0;
    int count = 0;
    
    if (_molecule->isPolymer())
    {
        PolymerPtr polymer = ToPolymerPtr(_molecule);
        
        for (int i = 0; i < polymer->monomerCount(); i++)
        {
            MonomerPtr monomer = polymer->getMonomer(i);
            
            if (!monomer)
            {
                continue;
            }

            ResButton *button = new ResButton(this, monomer);
            button->setGeometry(currX, 0, BUTTON_DIM, BUTTON_DIM);
            button->show();
	    connect(button, SIGNAL(clicked()), this,
                SLOT(pushResidueButton()));
            currX += BUTTON_DIM;
	    _buttons.push_back(button);
        }
    }
    
    setMinimumSize(QSize(currX, 30));
    setGeometry(0, 0, currX, 30);
}

void SequenceView::pushResidueButton()
{
    QObject *object = sender();
    ResButton *button = static_cast<ResButton *>(object);

    MonomerPtr monomer = button->getMonomer();
 
    _explorer->displayMonomer(monomer);

}

SequenceView::~SequenceView()
{
    for (int i = 0; i < _buttons.size(); i++)
    {
        delete _buttons[i];
        _buttons[i] = 0;
    }
    
    _buttons.clear();
}
