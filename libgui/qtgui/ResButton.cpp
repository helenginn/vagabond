#include "ResButton.h"
#include "../../libsrc/Monomer.h"

void ResButton::getText()
{   
    std::string code = _monomer->getResCode();
    QString qCode = QString::fromStdString(code);
    setText(qCode);
}
