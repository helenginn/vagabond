#include "Group.h"
#include <fstream>
#include "MtzFFT.h"
#include "MtzFile.h"
#include <libsrc/Crystal.h>
#include <libsrc/CSV.h>
#include <libsrc/Atom.h>
#include <libsrc/Monomer.h>
#include <libsrc/Polymer.h>
#include <libica/svdcmp.h>
#include "FileReader.h"
#include "MatrixView.h"
#include "AveDiffraction.h"
#include "AveCAlpha.h"
#include "AveUnitCell.h"

#include <QObject>
#include <QColor>

Group *Group::_topGroup = NULL;

Group::Group(QTreeWidget *parent) : QTreeWidgetItem(parent)
{
	_which = GroupMe;
	_exported = false;
	_centre = empty_vec3();
	_type = AveDiff;
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
	_mySet.recip = NULL;
	_mySet.ca = NULL;
	_mySet.unitCell = NULL;
	
	if (_topGroup == NULL)
	{
		_topGroup = this;
	}
}

void Group::performAverage()
{
	_message = "";

	std::cout << "Calculating averages." << std::endl;
	calculateAllAverages(true);
	std::cout << "Average sorted" << std::endl;
	
	if (_message.length() > 0)
	{
		emit failed();
	}
	else
	{
		emit resultReady();
	}
}

void Group::calculateAllAverages(bool force)
{
	if (force || _mySet.recip == NULL)
	{
		delete _mySet.recip;
		_mySet.recip = new AveDiffraction(this, _res);
		_mySet.recip->calculate();
	}

	if (force || _mySet.ca == NULL)
	{
		delete _mySet.ca;
		_mySet.ca = new AveCAlpha(this);
		_mySet.ca->calculate();
	}

	if (force || _mySet.unitCell == NULL)
	{
		delete _mySet.unitCell;
		_mySet.unitCell = new AveUnitCell(this);
		_mySet.unitCell->calculate();
	}
}

void Group::performCluster()
{
	_message = "";
	std::cout << "Calculating averages from " << mtzCount() <<
	" datasets..." << std::endl;
	calculateAllAverages(false);
	std::cout << "Inter-correlations against average..." << std::endl;
	findIntercorrelations();
	std::cout << "Calculating SVD..." << std::endl;
	svd();
	drawAxes();

	if (_message.length() > 0)
	{
		std::cout << "Failed clustering. " << std::endl;
		emit failed();
	}
	else
	{
		std::cout << "Clustering success." << std::endl;
		emit resultReady();
	}
}

void Group::copyFromOriginal(Group *ave)
{
	std::cout << "Copying from original..." << std::endl;
	_origSet = *(ave->getWorkingSet());
	_type = ave->_type;
}

void Group::updateText()
{
	std::string str = "Group of " + i_to_str(_mtzs.size()) + " data sets";
	setText(0, QString::fromStdString(str));
}

void Group::addMtz(MtzFFTPtr mtz)
{
	MtzFFTPtr tmp = MtzFFTPtr(new MtzFFT(this, *mtz));
	updateText();
	_mtzs.push_back(tmp);
}

void Group::addMtz(DiffractionMtzPtr mtzDiff, MtzFile *file,
                      CrystalPtr crystal)
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
	tmp->setSpaceGroup(mtzDiff->getFFT()->getSpaceGroup());
	tmp->setMtzFile(file);
	tmp->wipe();
	tmp->updateText();
	
	file->setCrystal(crystal);
	
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

AverageSet *Group::getWorkingSet()
{
	AverageSet *working = &_mySet;
	
	if (_which == GroupTop)
	{
		std::cout << "Using top group average" << std::endl;
		working = &_topGroup->_mySet;
	}
	else if (_which == GroupOriginal)
	{
		std::cout << "Using original group average" << std::endl;
		working = &_origSet;
	}
	else
	{
		std::cout << "Using my group average" << std::endl;
	}

	return working;
}

void Group::findIntercorrelations()
{
	svdAlloc();
	
	AverageSet *working = getWorkingSet();

	if (_type == AveDiff)
	{
		working->recip->findIntercorrelations(this, _svdPtrs);
	}
	else if (_type == AveCA)
	{
		working->ca->findIntercorrelations(this, _svdPtrs);
	}
	else if (_type == AveUC)
	{
		working->unitCell->findIntercorrelations(this, _svdPtrs);
	}

	size_t dims = _mtzs.size();
	memcpy(_orig, _svd, sizeof(double) * dims * dims);

	_correlMatrix = new MatrixView(this, dims, dims);
	_correlMatrix->populate();
}

void Group::drawResults(double **data, std::string filename)
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

void Group::svdAlloc()
{
	cleanupSVD();

	matAlloc(&_svd, &_svdPtrs);
	matAlloc(&_orig, &_origPtrs);
	matAlloc(&_v, &_vPtrs);

	size_t svd_dims = _mtzs.size();
	_w = (double *)malloc(sizeof(double) * svd_dims);
	memset(_w, '\0', sizeof(double) * svd_dims);
}

void Group::matAlloc(double **raw, double ***ptrs)
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

void Group::svd()
{
	size_t dims = _mtzs.size();

	int success = svdcmp((mat)_svdPtrs, dims, dims, (vect)_w, (mat)_vPtrs);

	if (success == 0)
	{
		_message = "Single value decomposition failed.";
		return;
	}
}

vec3 Group::getPoint(int num, int a1, int a2, int a3)
{
	vec3 v;
	v.x = _clusterPtrs[num][a1];
	v.y = _clusterPtrs[num][a2];
	v.z = _clusterPtrs[num][a3];

	return v;
}

void Group::drawAxes()
{
	matAlloc(&_cluster, &_clusterPtrs);

	for (size_t i = 0; i < _mtzs.size(); i++)
	{
		for (size_t j = 0; j < _mtzs.size(); j++)
		{
			_clusterPtrs[i][j] = _w[j] * (_svdPtrs[i][j]);
		}
	}
}

void Group::cleanupSVD()
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

void Group::setMarked(bool marked)
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

void Group::setDead(bool dead)
{
	if (_topGroup == this)
	{
		return;
	}

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

void Group::writeToStream(std::ofstream &f, bool complete)
{
	if (complete && !_marked)
	{
		return;
	}

	for (size_t i = 0; i < _mtzs.size(); i++)
	{
		MtzFile *file = _mtzs[i]->getMtzFile();

		f << file->getFilename();

		if (i < _mtzs.size() - 1)
		{
			f << ",";
		}
	}

	f << std::endl;
}

MtzFile *Group::getMtzFile(int i)
{
	return _mtzs[i]->getMtzFile();
}

void Group::flipMtzSelection(int i)
{
	MtzFile *file = _mtzs[i]->getMtzFile();
	if (!file) 
	{
		return;
	}

	file->flipSelected();
}

void Group::setMtzSelection(size_t i, bool val)
{
	if (i >= _mtzs.size())
	{
		return;
	}

	MtzFile *file = _mtzs[i]->getMtzFile();
	file->setSelected(val);
}

VagFFTPtr Group::getAverageFFT()
{
	return _mySet.recip->getFFT();
}

vec3 Group::getCentre()
{
	return _topGroup->_mySet.ca->getCentre();
}

Group::~Group()
{
	if (this == _topGroup)
	{
		_topGroup = NULL;
	}
}

void Group::averageRs(double *rwork, double *rfree, 
                      double *swork, double *sfree)
{
	CorrelData cd = empty_CD();

	for (size_t i = 0; i < mtzCount(); i++)
	{
		MtzFile *file = getMtzFile(i);
		
		if (file->isDead())
		{
			continue;
		}
		
		add_to_CD(&cd, file->getRWork(), file->getRFree());
	}
	
	means_stdevs_CD(cd, rwork, rfree, swork, sfree);
}

std::vector<MtzFFTPtr> Group::getMtzsFromSelection()
{
	std::vector<MtzFFTPtr> mtzs;
	for (size_t i = 0; i < mtzCount(); i++)
	{
		MtzFile *file = getMtz(i)->getMtzFile();
		if (file->isSelected())
		{
			mtzs.push_back(getMtz(i));
		}
	}
	
	return mtzs;
}

std::vector<double> Group::getUnitCell()
{
	return getWorkingSet()->unitCell->getUnitCell();
}

void Group::changeColour(double r, double b, double g, double a)
{
	for (size_t i = 0; i < mtzCount(); i++)
	{
		MtzFile *file = getMtz(i)->getMtzFile();
		file->setColour(r, g, b, a);
	}
}
