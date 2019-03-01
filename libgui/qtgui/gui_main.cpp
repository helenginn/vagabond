//
//  main.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "../../libsrc/Options.h"
#include <execinfo.h>
#include <signal.h>
#include "StartScreen.h"
#include <QtWidgets/qapplication.h>

void handler(int sig) {
	void *array[10];
	size_t size;

	// get void*'s for all entries on the stack
	size = backtrace(array, 10);

	// print out all the frames to stderr
	fprintf(stderr, "Error: signal %d:\n", sig);
	backtrace_symbols_fd(array, (int)size, STDERR_FILENO);
	exit(1);
}

int main(int argc, char * argv[])
{
	// insert code here...

	signal(SIGSEGV, handler);
	signal(SIGABRT, handler);

	std::cout << "Qt version: " << qVersion() << std::endl;

	QApplication app(argc, argv);
	setlocale(LC_NUMERIC, "C");

	StartScreen startScreen(NULL, argc, argv);
	startScreen.show();
	int status = 0;

	status = app.exec();
	
	if (status == 0)
	{
//		InstructionThread *thread = window.getInstructionThread();
	}
	
	return status;
}
