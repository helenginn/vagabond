#include "InstructionThread.h"
#include "VagWindow.h"

void InstructionThread::run()
{
    _window->waitForInstructions();
}