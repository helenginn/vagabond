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
#include "Molecule.h"
#include <math.h>
#include <iomanip>

int main(int argc, const char * argv[]) {
	// insert code here...
	std::cout << "Hello, World!\n";

	MoleculePtr mol;

	PDBReader pdb = PDBReader();
	pdb.setFilename("5i40_final.pdb");
	mol = pdb.getMolecule();

	FFTPtr fft = FFTPtr(new cFFTW3d());
	mol->calculateMillers(fft);

	fft->createFFTWplan(8);

	fft->fft(1);

	for (int i = -24; i < 24; i++)
	for (int j = -34; j < 34; j++)
	for (int k = -40; k < 40; k++)
	{
		std::cout << std::fixed << std::setprecision(1) << std::setw(4) << i << std::setw(4) << j << std::setw(4) << k << std::setw(8) << std::right << sqrt(fft->getIntensity(i, j, k)) <<  " 1.0000  " << std::setw(5) << std::right << fft->getPhase(i, j, k) << std::setw(8) << 1000 << std::endl;
	}


    return 0;
}
