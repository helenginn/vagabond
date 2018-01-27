#ifndef __Vagabond__InstructionThread__
#define __Vagabond__InstructionThread__

#include <QtCore/qthread.h>

class VagWindow;

class InstructionThread : public QThread
{
public:
    InstructionThread() : QThread(NULL)
    {
    
    }

	InstructionThread(QObject *object) : QThread(object)
	{
	
	}
    
    void run();
    
    void setVagWindow(VagWindow *window)
    {
        _window = window;
    }

private:
    VagWindow *_window;
};

#endif