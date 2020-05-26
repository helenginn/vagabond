#include <iostream>
#include <iomanip>
#include <libsrc/shared_ptrs.h>
#include <libsrc/DiffractionMTZ.h>
#include <libsrc/Shouter.h>
#include <libsrc/FileReader.h>
#include <libsrc/Crystal.h>
#include <libsrc/FFT.h>

int main(int argc, char * argv[])
{
	std::cout << "Shell-by-shell scaling." << std::endl;

	if (argc <= 2)
	{
		std::cout << "Not enough arguments, please supply MTZ files!" 
		<< std::endl << std::endl;

		std::cout << "Minimum form required to execute:" << std::endl;
		std::cout << "shell_scale reference.mtz change_me.mtz" << std::endl;
		std::cout << std::endl;
	}

	DiffractionMtzPtr ref_ = DiffractionMtzPtr(new DiffractionMtz());
	DiffractionMtzPtr mtz_ = DiffractionMtzPtr(new DiffractionMtz());

	try
	{
		ref_->setNeedsRfree(false);
		ref_->setFilename(std::string(argv[1]));
		ref_->load();

		mtz_->setNeedsRfree(false);
		mtz_->setFilename(std::string(argv[2]));
		mtz_->load();
	}
	catch (Shouter *s)
	{
		s->shoutToStdOut();
	}

	VagFFTPtr ref = ref_->getFFT();
	VagFFTPtr mtz = mtz_->getFFT();
	mtz->setStatus(FFTReciprocalSpace);
	
	double maxref = ref_->getMaxResolution();
	double maxmtz = mtz_->getMaxResolution();
	
	if (maxref > maxmtz)
	{
		std::cout << "WARNING!! Resolution of reference (" << maxref
		<< " Å) is worse than that of the dataset to be scaled "
		<< "( " << maxmtz << " Å)." << std::endl << " This means that "
		<< " some reflections cannot be scaled against the reference."
		<< std::endl << "To stop things going wrong, unmatched "
		<< " reflections will be set to zero." << std::endl;
	}


	std::string whole = std::string(argv[2]);
	std::string path = getPath(whole);
	std::string file = getFilename(whole);
	std::string newfile = path + "shsc-" + file;

	std::vector<ShellInfo> _shells;
	makeShells(&_shells, 0, maxref);

	vec3 nLimits = getNLimits(ref, mtz);

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				double res = ref->resolution(i, j, k);

				int s = findShell(_shells, res);
				
				if (s < 0)
				{
					continue;
				}

				double fr = ref->getReal(i, j, k);
				double fm = mtz->getReal(i, j, k);

				_shells[s].work1.push_back(fr);
				_shells[s].work2.push_back(fm);
			}
		}
	}

	
	for (size_t i = 0; i < _shells.size(); i++)
	{
		_shells[i].scale = scale_factor_by_sum(_shells[i].work1,
		                                       _shells[i].work2);
		_shells[i].rFactor = r_factor(_shells[i].work1, _shells[i].work2);
	}
	
	std::cout << std::endl;

	std::cout << std::setprecision(3);
	std::cout << std::setw(19) << "Resolution bin  " << std::flush;
	std::cout << std::setw(13) << "R factor" << std::flush;
	std::cout << std::setw(15) << "Scale factor" << std::flush;
	std::cout << std::endl;

	for (size_t i = 0; i < _shells.size(); i++)
	{
		ShellInfo *sh = &_shells[i];
		std::cout << std::setw(6) << sh->maxRes << " - " << 
		std::setw(6) << sh->minRes << " Å";
		std::cout << std::setw(15) << sh->rFactor;
		std::cout << std::setw(15) << sh->scale;
		std::cout << std::endl;
	}
	
	std::cout << std::endl;
	
	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				double res = mtz->resolution(i, j, k);
				int s = findShell(_shells, res);
				
				if (s < 0)
				{
					mtz->setReal(i, j, k, 0);
					continue;
				}

				double scale = _shells[s].scale;
				double amp = mtz->getReal(i, j, k);
				double sigma = mtz->getReal(i, j, k);
				amp /= scale;
				sigma /= scale;
				mtz->setReal(i, j, k, amp);
				mtz->setImag(i, j, k, sigma);
			}
		}
	}
	
	mtz->writeToFile(newfile, maxmtz, mtz);
	std::cout << "Written MTZ file to " << newfile << "." << std::endl;
	

}
