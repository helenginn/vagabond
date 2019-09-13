#include "InstructionThread.h"
#include "VagWindow.h"

void InstructionThread::run()
{
    int exitCode = _window->waitForInstructions();

}
