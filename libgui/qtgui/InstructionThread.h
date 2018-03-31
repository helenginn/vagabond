#ifndef __Vagabond__InstructionThread__
#define __Vagabond__InstructionThread__

#include <QtCore/qthread.h>

class VagWindow;

class InstructionThread : public QThread
{
public:
    InstructionThread() : QThread(NULL)
    {
    	_mustDie = false;
    }

	InstructionThread(QObject *object) : QThread(object)
	{
    	_mustDie = false;
	}
    
    void run();
    
    void setVagWindow(VagWindow *window)
    {
        _window = window;
    }

	void setShouldDie()
	{
		_mustDie = true;
	}
	
	bool shouldDie()
	{
		return _mustDie;	
	}
private:
    VagWindow *_window;
	bool _mustDie;
};

#endif
