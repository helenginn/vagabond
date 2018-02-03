#include "../../libsrc/Options.h"
#include "../../libsrc/Notifiable.h"
#include "SetterEdit.h"
#include <iostream>

void SetterEdit::setObjectValue()
{
    QString myText = text();
    std::string stdText = myText.toStdString();
    float value = atof(stdText.c_str());

    if (_degrees)
    {
        value = deg2rad(value);
    }

    OptionsPtr options = Options::getRuntimeOptions();
    Notifiable *notify = options->getNotify();

    std::cout << "Got here." << std::endl;

    notify->setObject(_object);
    notify->setSetter(_setter, value);
    notify->setRefreshGroup(_monomer);
    notify->setInstruction(InstructionTypeSetObjectValue);
}


