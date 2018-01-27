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

    (*_setter)(_object, value);

    if (_monomer)
    {   
        _monomer->propagateChange();
        _monomer->refreshPositions();
    }
}


