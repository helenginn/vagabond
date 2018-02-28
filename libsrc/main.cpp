//
//  main.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Options.h"
#include <execinfo.h>
#include <signal.h>


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



int main(int argc, const char * argv[])
{
	/* Options are parsed and will generate the objects needed */

	signal(SIGSEGV, handler);
	signal(SIGABRT, handler);


	Options options(argc, argv);
	options.run();

	return 0;
}
