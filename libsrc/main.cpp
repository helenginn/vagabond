//
//  main.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Options.h"

int main(int argc, const char * argv[])
{
	/* Options are parsed and will generate the objects needed */

	Options options(argc, argv);
	options.run();

    return 0;
}
