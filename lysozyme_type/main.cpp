#include <iostream>
#include <fstream>
#include <libsrc/shared_ptrs.h>
#include <libsrc/DiffractionMTZ.h>
#include <libsrc/FFT.h>
#include <libsrc/Shouter.h>
#include <hcsrc/FileReader.h>

int main(int argc, char * argv[])
{
	std::cout << "Forcing the Wilson plot down." << std::endl;

	if (argc <= 1)
	{
		std::cout << "Not enough arguments, please supply MTZ files!" << std::endl;
	}
	
	std::ofstream file;
	file.open("lysozyme_types.csv");
	file << "filename, 11_11_4, 11_11_5, ratio" << std::endl;
	for (size_t i = 1; i < argc; i++)
	{
		DiffractionMtzPtr mtz = DiffractionMtzPtr(new DiffractionMtz());

		try
		{
			mtz->setNeedsRfree(false);
			mtz->setFilename(std::string(argv[i]));
			mtz->load();
			
			VagFFTPtr fft = mtz->getFFT();
			float f11_11_4 = fft->getReal(11, 11, 4);
			float f11_11_5 = fft->getReal(11, 11, 5);
			
			file << argv[i] << ", " << f11_11_4 << ", " << f11_11_5 << ", " << f11_11_4 / f11_11_5 << std::endl;


		}
		catch (Shouter *s)
		{
			s->shoutToStdOut();
		}

	}
	
	file.close();
	
	std::cout << "Written to lysozyme_types.csv" << std::endl;
}
