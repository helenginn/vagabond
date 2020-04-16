#include "Averager.h"
#include <fstream>
#include "MtzFFT.h"
#include "MtzFile.h"
#include <libsrc/Crystal.h>
#include <libsrc/CSV.h>
#include <libica/svdcmp.h>
#include "FileReader.h"
#include "MatrixView.h"

Averager::Averager(QTreeWidget *parent) : QTreeWidgetItem(parent)
{
	_res = 2;
	_svd = NULL;
	_svdPtrs = NULL;
	_orig = NULL;
	_origPtrs = NULL;
	_v = NULL;
	_vPtrs = NULL;
	_w = NULL;
	_correlMatrix = NULL;
	_marked = false;
}

void Averager::performAverage()
{
	_message = "";

	makeAverage(true);
	
	if (_message.length() > 0)
	{
		emit failed();
	}
	else
	{
		emit resultReady();
	}
}

void Averager::performCluster()
{
	_message = "";
	std::cout << "Calculating average from " << mtzCount() <<
	" datasets..." << std::endl;
	makeAverage();
	std::cout << "Scaling individual MTZs to average..." << std::endl;
	scaleIndividuals();
	std::cout << "Inter-correlations against average..." << std::endl;
	findIntercorrelations();
	std::cout << "Calculating SVD..." << std::endl;
	svd();
	drawAxes();

	if (_message.length() > 0)
	{
		emit failed();
	}
	else
	{
		emit resultReady();
	}
}

void Averager::useOriginalAverage()
{
	if (_origAve)
	{
		_fft = _origAve;
	}
}

void Averager::copyFromAverage(Averager *ave)
{
	_fft = ave->_fft;
}

void Averager::makeAverage(bool force)
{
	if (_mtzs.size() == 0)
	{
		_message += "No data sets in this group.\n";
		return;
	}
	
	if (_fft && _fft->nn() > 0 && !force)
	{
		return;
	}
	
	if (force)
	{
		_origAve = VagFFTPtr(new VagFFT(*_fft));
	}

	VagFFTPtr first = _mtzs[0];
	_fft = VagFFTPtr(new VagFFT(*first));
	_fft->wipe();
	int count = 0;

	for (size_t i = 0; i < _mtzs.size(); i++)
	{
		MtzFFTPtr current = _mtzs[i];
		if (current->getMtzFile()->isDead())
		{
			continue;
		}
		
		count++;

		vec3 nLimits = getNLimits(_fft, current);

		for (int k = -nLimits.z; k < nLimits.z; k++)
		{
			for (int j = -nLimits.y; j < nLimits.y; j++)
			{
				for (int i = -nLimits.x; i < nLimits.x; i++)
				{
					double res = current->resolution(i, j, k);
					if (res < _res)
					{
						continue;
					}

					double amp = current->getReal(i, j, k);
					long ele = _fft->element(i, j, k);
					if (amp == amp)
					{
						_fft->addToReal(ele, amp);
					}
				}
			}
		}
	}
	
	_fft->multiplyAll(1 / (double)count);
}

void Averager::updateText()
{
	std::string str = "Group of " + i_to_str(_mtzs.size()) + " data sets";
	setText(0, QString::fromStdString(str));
}

void Averager::addMtz(MtzFFTPtr mtz)
{
	MtzFFTPtr tmp = MtzFFTPtr(new MtzFFT(this, *mtz));
	updateText();
	_mtzs.push_back(tmp);
}

void Averager::addMtz(DiffractionMtzPtr mtzDiff, MtzFile *file)
{
//	_names.push_back(mtzDiff->getFilename());

	VagFFTPtr ref = mtzDiff->getFFT();

	if (_mtzs.size() > 0)
	{
		ref = _mtzs[0];
	}

	/* ensure identical spacing for subsequent mtzs added */
	MtzFFTPtr tmp = MtzFFTPtr(new MtzFFT(this, *ref));
	std::vector<double> uc = mtzDiff->getFFT()->getUnitCell();
	tmp->setUnitCell(uc);
	tmp->setText(0, QString::fromStdString(mtzDiff->getFilename()));
	tmp->setMtzFile(file);
	tmp->wipe();
	
	VagFFTPtr mtz = mtzDiff->getFFT();

	vec3 nLimits = getNLimits(tmp, mtz);

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				double amp = mtz->getReal(i, j, k);
				long ele = tmp->element(i, j, k);
				tmp->addToReal(ele, amp);
			}
		}
	}
	
	_mtzs.push_back(tmp);
	updateText();
}

void Averager::scaleIndividualMtz(int i)
{
	VagFFTPtr one = _mtzs[i];
	std::vector<ShellInfo> _shells;
	makeShells(&_shells, 0, _res);

	vec3 nLimits = getNLimits(_fft, one);

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				double res = one->resolution(i, j, k);
				if (res < _res)
				{
					continue;
				}

				double amp = one->getReal(i, j, k);
				double ref = _fft->getReal(i, j, k);

				int s = findShell(_shells, res);
				_shells[s].work1.push_back(amp);
				_shells[s].work2.push_back(ref);
			}
		}
	}
	
	for (size_t i = 0; i < _shells.size(); i++)
	{
		_shells[i].scale = scale_factor_by_sum(_shells[i].work1,
		                                       _shells[i].work2);
	}
	
	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				double res = one->resolution(i, j, k);
				int s = findShell(_shells, res);
				
				if (s < 0)
				{
					continue;
				}

				double scale = _shells[s].scale;
				double amp = one->getReal(i, j, k);
				amp *= scale;
				one->setReal(i, j, k, amp);
			}
		}
	}
}

void Averager::scaleIndividuals()
{
	for (size_t i = 0; i < _mtzs.size(); i++)
	{
		scaleIndividualMtz(i);
	}
}

void Averager::findIntercorrelations()
{
	svdAlloc();

	for (size_t i = 1; i < _mtzs.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			double cc = 0;
			
			if (i != j)
			{
				cc = findCorrelation(i, j);
			}
			
			_svdPtrs[i][j] = cc;
			_svdPtrs[j][i] = cc;
		}

		std::cout << "." << std::flush;
	}

	std::cout << std::endl;

	size_t dims = _mtzs.size();
	memcpy(_orig, _svd, sizeof(double) * dims * dims);

	_correlMatrix = new MatrixView(this, dims, dims);
	_correlMatrix->populate();
	drawResults(_svdPtrs, "correlation_0");
}

void Averager::drawResults(double **data, std::string filename)
{

	CSVPtr csv = CSVPtr(new CSV(3, "i", "j", "cc"));
	
	for (size_t i = 0; i < _mtzs.size(); i++)
	{
		for (size_t j = 0; j < _mtzs.size(); j++)
		{
			double cc = data[i][j];
			csv->addEntry(3, (double)i, (double)j, (double)cc);
		}
	}

	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = filename;
	plotMap["height"] = "800";
	plotMap["width"] = "800";
	plotMap["xHeader0"] = "i";
	plotMap["yHeader0"] = "j";
	plotMap["zHeader0"] = "cc";

	plotMap["xTitle0"] = "mtz num";
	plotMap["yTitle0"] = "mtz num";
	plotMap["style0"] = "heatmap";
	plotMap["stride"] = i_to_str(_mtzs.size());

	csv->plotPNG(plotMap);

}

VagFFTPtr Averager::makeDifference(VagFFTPtr one, VagFFTPtr two)
{
	std::vector<double> xs, ys;
	VagFFTPtr fft = VagFFTPtr(new VagFFT(*one));

	for (int i = 0; i < one->nn(); i++)
	{
		double yours = two->getReal(i);

		fft->addToReal(i, -yours);
	}

	return fft;
}

double Averager::pairwiseValue(VagFFTPtr oneTwoDiff,
                               VagFFTPtr three, VagFFTPtr four)
{
	CorrelData cd = empty_CD();

	for (int i = 0; i < oneTwoDiff->nn(); i++)
	{
		double x = oneTwoDiff->getReal(i);
		double y = three->getReal(i);
		y -= four->getReal(i);

		add_to_CD(&cd, x, y);
	}

	double cc = evaluate_CD(cd);
	if (cc < 0) cc = 0;
	return cc;
}

void Averager::exhaustiveCorrelations()
{
	for (size_t i = 1; i < _mtzs.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			VagFFTPtr one = _mtzs[i];
			VagFFTPtr two = _mtzs[j];
			VagFFTPtr oneTwoDiff = makeDifference(one, two);
			
			Pair oneTwo = makePair(one, two);
			Pair twoOne = makePair(two, one);

			for (size_t k = 1; k < j; k++)
			{
				for (size_t l = 0; l < k; l++)
				{
					VagFFTPtr three = _mtzs[k];
					VagFFTPtr four = _mtzs[l];

					Pair threeFour = makePair(three, four);
					Pair fourThree = makePair(four, three);

					if (_vagVals.count(oneTwo) == 0 || _vagVals.count(twoOne) == 0)
					{
						_vagVals[oneTwo][three] = 0;
						_vagVals[twoOne][three] = 0;
						_vagVals[oneTwo][four] = 0;
						_vagVals[twoOne][four] = 0;
					}

					if (_vagVals.count(threeFour) == 0 || 
					    _vagVals.count(fourThree) == 0)
					{
						_vagVals[threeFour][one] = 0;
						_vagVals[fourThree][one] = 0;
						_vagVals[threeFour][two] = 0;
						_vagVals[fourThree][two] = 0;
					}

					double diff = pairwiseValue(oneTwoDiff, three, four);
					
					_vagVals[oneTwo][three] += diff;
					_vagVals[twoOne][three] += diff;
					_vagVals[oneTwo][four] += diff;
					_vagVals[twoOne][four] += diff;
					
					_vagVals[threeFour][one] += diff;
					_vagVals[threeFour][one] += diff;
					_vagVals[fourThree][two] += diff;
					_vagVals[fourThree][two] += diff;
				}
			}
		}
		
		std::cout << "." << std::flush;
	}
	std::cout << std::endl;
	
	std::cout << "Pairwise comparisons" << std::endl;
	svdAlloc();

	for (size_t i = 1; i < _mtzs.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			double cc = 0;
			
			if (i != j)
			{
				cc = comparePairwise(i, j);
			}
			
			_svdPtrs[i][j] = cc;
			_svdPtrs[j][i] = cc;
		}
		std::cout << "." << std::flush;
	}

	std::cout << std::endl;
	drawResults(_svdPtrs, "correlation_0");
}

double Averager::comparePairwise(int i, int j)
{
	VagFFTPtr one = _mtzs[i];
	VagFFTPtr two = _mtzs[j];

	CorrelData cd = empty_CD();
		
	for (size_t i = 0; i < _mtzs.size(); i++)
	{
		{
			Pair withOne = makePair(one, _mtzs[i]);
			Pair withTwo = makePair(two, _mtzs[i]);

			VagScore::iterator it;
			for (it = _vagVals[withOne].begin(); 
			     it != _vagVals[withOne].end(); it++)
			{
				VagFFTPtr pair = it->first;
				double c1 = it->second;

				if (_vagVals[withTwo].count(pair) == 0)
				{
					continue;
				}

				double c2 = _vagVals[withTwo][it->first];
				add_to_CD(&cd, c1, c2);

			}
		}
	}


	double cc = evaluate_CD(cd);
	return cc;
}

double Averager::findCorrelation(int i, int j)
{
	VagFFTPtr one = _mtzs[i];
	VagFFTPtr two = _mtzs[j];

	CorrelData cd = empty_CD();
	
	for (int i = 0; i < one->nn(); i++)
	{
		double ave = _fft->getReal(i);
		double mine = one->getReal(i);
		double yours = two->getReal(i);

		double x = mine - ave;
		double y = yours - ave;

		add_to_CD(&cd, x, y);
	}
	
	double cc = evaluate_CD(cd);
	if (cc < 0) return 0;
	return cc;
}

void Averager::svdAlloc()
{
	cleanupSVD();

	matAlloc(&_svd, &_svdPtrs);
	matAlloc(&_orig, &_origPtrs);
	matAlloc(&_v, &_vPtrs);
	
	size_t svd_dims = _mtzs.size();
	_w = (double *)malloc(sizeof(double) * svd_dims);
	memset(_w, '\0', sizeof(double) * svd_dims);
}

void Averager::matAlloc(double **raw, double ***ptrs)
{
	size_t svd_dims = _mtzs.size();
	size_t square = svd_dims * svd_dims;

	*raw = (double *)malloc(sizeof(double *) * square);
	memset(*raw, '\0', sizeof(double) * square);

	*ptrs = (double **)malloc(sizeof(double **) * svd_dims);

	for (size_t i = 0; i < svd_dims; i++)
	{
		(*ptrs)[i] = &((*raw)[i * svd_dims]);
	}
}

void Averager::svd()
{
	size_t dims = _mtzs.size();

	int success = svdcmp((mat)_svdPtrs, dims, dims, (vect)_w, (mat)_vPtrs);
	
	if (success == 0)
	{
		_message = "Single value decomposition failed.";
		return;
	}
	
	drawResults(_svdPtrs, "correlation_1");
}

vec3 Averager::getPoint(int num, int a1, int a2, int a3)
{
	vec3 v;
	v.x = _clusterPtrs[num][a1];
	v.y = _clusterPtrs[num][a2];
	v.z = _clusterPtrs[num][a3];

	return v;
}

void Averager::drawAxes()
{
	CSVPtr csv = CSVPtr(new CSV(7, "mtz", "a", "b", "c", "d", "e", "f"));

	matAlloc(&_cluster, &_clusterPtrs);

	for (size_t i = 0; i < _mtzs.size(); i++)
	{
		double *vec = _origPtrs[i];
		double *aligned = (double *)malloc(sizeof(double) * _mtzs.size());
		memset(aligned, '\0', sizeof(double) * _mtzs.size());

		for (size_t j = 0; j < _mtzs.size(); j++)
		{
			aligned[j] += vec[j] * _svdPtrs[i][j];
			aligned[j] *= _w[j];
			_clusterPtrs[i][j] = aligned[j];
		}
		
		double dist = 0;
		for (int i = 0; i < 5; i++)
		{
			dist += aligned[0] * aligned[0];
		}
		
		dist = sqrt(dist);
		
		csv->addEntry(7, (double)i, aligned[0], aligned[1], aligned[2],
		              aligned[3], aligned[4], aligned[5], aligned[6]);
		
		free(aligned);
	}
	
	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = "principal_axes";
	plotMap["xHeader0"] = "a";
	plotMap["yHeader0"] = "b";
	plotMap["colour0"] = "black";

	plotMap["xTitle0"] = "p-ness";
	plotMap["yTitle0"] = "q-ness";
	plotMap["style0"] = "scatter";
	csv->plotPNG(plotMap);

	csv->writeToFile("principal_axes.csv");
}

void Averager::cleanupSVD()
{
	free(_svd);
	free(_svdPtrs);
	free(_orig);
	free(_origPtrs);
	free(_v);
	free(_vPtrs);
	free(_w);
	_svd = NULL;
	_svdPtrs = NULL;
	_orig = NULL;
	_origPtrs = NULL;
	_cluster = NULL;
	_clusterPtrs = NULL;
	_v = NULL;
	_vPtrs = NULL;
	_w = NULL;
}

void Averager::setMarked(bool marked)
{
	for (size_t i = 0; i < _mtzs.size(); i++)
	{
		MtzFile *file = _mtzs[i]->getMtzFile();
		file->setMarked(marked);
	}

	_marked = marked;

	QFont curr = font(0);
	curr.setBold(marked);
	setFont(0, curr);
}

void Averager::setDead(bool dead)
{
	for (size_t i = 0; i < _mtzs.size(); i++)
	{
		MtzFile *file = _mtzs[i]->getMtzFile();
		file->setDead(dead);
	}

	_dead = dead;

	QFont curr = font(0);
	curr.setItalic(_dead);
	setFont(0, curr);
}

void Averager::writeToStream(std::ofstream &f, bool complete)
{
	for (size_t i = 0; i < _mtzs.size(); i++)
	{
		MtzFile *file = _mtzs[i]->getMtzFile();
		
		if (complete && !file->isMarked())
		{
			continue;
		}
		
		f << file->getFilename();
		
		if (i < _mtzs.size() - 1)
		{
			f << ",";
		}
	}

	f << std::endl;
}

MtzFile *Averager::getMtzFile(int i)
{
	return _mtzs[i]->getMtzFile();
}

void Averager::flipMtzSelection(int i)
{
	MtzFile *file = _mtzs[i]->getMtzFile();
	if (!file) 
	{
		return;
	}
	
	file->flipSelected();
}

void Averager::setMtzSelection(size_t i, bool val)
{
	if (i >= _mtzs.size())
	{
		return;
	}

	MtzFile *file = _mtzs[i]->getMtzFile();
	file->setSelected(val);
}