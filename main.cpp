//
//  main.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include <iostream>
#include "libsrc/fftw3d.h"
#include "PDBReader.h"
#include "shared_ptrs.h"
#include "Crystal.h"
#include <math.h>
#include <iomanip>

int main(int argc, const char * argv[]) {
	// insert code here...
	std::cout << "Hello, World!\n";

	CrystalPtr crystal;

	PDBReader pdb = PDBReader();
	pdb.setFilename("5i40_final.pdb");
	crystal = pdb.getCrystal();

	FFTPtr fft = FFTPtr(new cFFTW3d());
	crystal->calculateMillers(fft);
	crystal->writeCalcMillersToFile(fft, 1.0);


	fft->fft(1);

    return 0;
}
