#include <iostream>
#include <libsrc/shared_ptrs.h>
#include <libsrc/DiffractionMTZ.h>
#include <libsrc/Shouter.h>
#include <hcsrc/FileReader.h>
#include "Wilson.h"

int main(int argc, char * argv[])
{
	std::cout << "Forcing the Wilson plot down." << std::endl;

	if (argc <= 1)
	{
		std::cout << "Not enough arguments, please supply MTZ file!" << std::endl;
	}

	DiffractionMtzPtr mtz = DiffractionMtzPtr(new DiffractionMtz());

	try
	{
		mtz->setNeedsRfree(false);
		mtz->setFilename(std::string(argv[1]));
		mtz->load();
	}
	catch (Shouter *s)
	{
		s->shoutToStdOut();
	}

	VagFFTPtr fft = mtz->getFFT();
	Wilson w(fft);
	std::string whole = std::string(argv[1]);
	std::string path = getPath(whole);
	std::string file = getFilename(whole);
	w.setFilename(path + "fd-" + file);
	w.forceDown();
}
