//
//  main.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "libsrc/fftw3d.h"
#include "PDBReader.h"
#include "shared_ptrs.h"
#include "Crystal.h"
#include "Options.h"
#include "Sandbox.h"

int main(int argc, const char * argv[])
{
	/* Options are parsed and will generate the objects needed */

	//outputTorsionStuff();

	Options options(argc, argv);
	options.run();

    return 0;
}
