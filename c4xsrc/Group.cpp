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
#include "QuickAtoms.h"
#include "AveCAlpha.h"
#include "AveCSV.h"
#include "AveUnitCell.h"

#include <QObject>
#include <QFontDatabase>

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
	_mySet.csv = NULL;
	_mySet.ca = NULL;
	_mySet.unitCell = NULL;
	Qt::ItemFlags fl = flags();
	setFlags(fl | Qt::ItemIsEditable);
	
	if (_topGroup == NULL)
	{
		_topGroup = this;
	}
	
}

void Group::setData(int column, int role, const QVariant &value)
{
	if (role == Qt::EditRole)
	{
		_customName = value.toString().toStdString();
	}
	
	if (_customName.length() == 0 && (int)value.type() == (int)QMetaType::QString)
	{
		std::string str = generateText();
		QVariant val = QVariant(QString::fromStdString(str));
		QTreeWidgetItem::setData(column, role, val);
		return;
	}

	QTreeWidgetItem::setData(column, role, value);

	return;
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
	
	if (force || _mySet.csv == NULL)
	{
		if (AveCSV::usingCSV())
		{
			delete _mySet.csv;
			_mySet.csv = new AveCSV(this, "");
			_type = AveComma;
		}
	}
}

void Group::performCluster()
{
	_message = "";
	std::cout << "Calculating averages from " << mtzCount() <<
	" datasets..." << std::endl;
	calculateAllAverages(true);
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

std::string Group::generateText()
{
	std::string str;
	
	if (_customName.length())
	{
		str += _customName;
	}
	else
	{
		str += "Group of " + i_to_str(_mtzs.size()) + " data sets";
	}

	return str;
}

void Group::updateText()
{
	std::string str = generateText();
	setText(0, QString::fromStdString(str));

	QFont curr = font(0);
	curr.setBold(_marked);
	setFont(0, curr);
}

void Group::addMtz(MtzFFTPtr mtz)
{
	MtzFFTPtr tmp = MtzFFTPtr(new MtzFFT(this, *mtz));
	updateText();
	_mtzs.push_back(tmp);
}

void Group::addMtz(DiffractionMtzPtr mtzDiff, MtzFile *file)
{
	VagFFTPtr ref = mtzDiff->getOriginal();

	if (!ref)
	{
		MtzFFTPtr tmp = MtzFFTPtr(new MtzFFT(this));
		tmp->setMtzFile(file);
		tmp->updateText();
		std::vector<double> uc = std::vector<double>(6, 0.);
		tmp->setUnitCell(uc);
		_mtzs.push_back(tmp);
		updateText();
		return;
	}

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
	
	VagFFTPtr mtz = mtzDiff->getOriginal();

	vec3 nLimits = getNLimits(tmp, mtz);

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				double amp = mtz->getReal(i, j, k);
				double imag = mtz->getImag(i, j, k);
				long ele = tmp->element(i, j, k);
				tmp->addToReal(ele, amp);
				tmp->addToImag(ele, imag);
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
		working = &_topGroup->_mySet;
	}
	else if (_which == GroupOriginal)
	{
		working = &_origSet;
	}
	else
	{

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
	else if (_type == AveComma)
	{
		if (working->csv != NULL)
		{
			working->csv->findIntercorrelations(this, _svdPtrs);
		}
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

/* formally known as void Group::svd() */
void Group::svd()
{
	size_t dims = _mtzs.size();
	
	for (size_t i = 0; i < dims; i++)
	{
		for (size_t j = 0; j < dims; j++)
		{
			if (_svdPtrs[i][j] != _svdPtrs[i][j])
			{
				_svdPtrs[i][j] = 0;
			}
		}
	}

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
	
	updateText();
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

void Group::writeToStream(std::ofstream &f, ExportType type, bool complete)
{
	if (complete && !_marked)
	{
		return;
	}

	for (size_t i = 0; i < _mtzs.size(); i++)
	{
		MtzFile *file = _mtzs[i]->getMtzFile();

		if (type == ExportNickname)
		{
			f << file->metadata();
		}
		else if (type == ExportFilename)
		{
			std::string fn = file->getFilename();
			f << getBaseFilename(fn);
		}
		else /* full path */
		{
			std::string fn = file->getFilename();
			f << getBaseFilenameWithPath(fn);
		}

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

std::string Group::getMetadata(int i)
{
	return getMtzFile(i)->metadata();
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
	if (_topGroup == NULL || _topGroup->_mySet.ca == NULL)
	{
		return empty_vec3();
	}

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
	if (getWorkingSet() != NULL && getWorkingSet()->unitCell != NULL)
	{
		return getWorkingSet()->unitCell->getUnitCell();
	}
	else 
	{
		return std::vector<double>(6, 0.);
	}
}

void Group::changeColour(double r, double b, double g, double a)
{
	for (size_t i = 0; i < mtzCount(); i++)
	{
		MtzFile *file = getMtz(i)->getMtzFile();
		file->setColour(r, g, b, a);
	}
}

void Group::useAverageType(GroupType type)
{
	if (type == AveComma && _mySet.csv == NULL)
	{
		return;
	}
	_type = type;
}

void Group::collapseDatasets()
{
	vec3 centre = empty_vec3();
	double count = 0;
	
	for (size_t i = 0; i < mtzCount(); i++)
	{
		if (getMtzFile(i)->isSelected())
		{
			vec3 c = getMtzFile(i)->getQuickAtoms()->getCentre();
			vec3_add_to_vec3(&centre, c);
			count++;
		}
	}
	
	vec3_mult(&centre, 1 / count);

	for (size_t i = 0; i < mtzCount(); i++)
	{
		MtzFile *file = getMtzFile(i);
		QuickAtoms *atoms = file->getQuickAtoms();
		atoms->collapseOnTarget(centre);
	}
	
	_topGroup->_mySet.ca = new AveCAlpha(this);
	_topGroup->_mySet.ca->calculate();
}

std::vector<MtzFFTPtr> Group::orderByCoverage()
{
	/* clusterPtrs(d, c); d = dataset, c = cluster */
	std::vector<MtzFFTPtr> reordered;
	double percentage = 0;
	
	size_t svd_dims = _mtzs.size();
	double *goals = (double *)malloc(sizeof(double) * svd_dims);
	memset(goals, '\0', sizeof(double) * svd_dims);

	double total = 0;
	for (size_t i = 0; i < mtzCount(); i++)
	{
		double max = 0;
		for (size_t j = 0; j < mtzCount(); j++)
		{
			max = std::max(max, _clusterPtrs[j][i]);
		}

		goals[i] = max;
		total += max;
	}
	
	double *covered = (double *)malloc(sizeof(double) * svd_dims);
	memset(covered, '\0', sizeof(double) * svd_dims);
	double *wip = (double *)malloc(sizeof(double) * svd_dims);
	memset(wip, '\0', sizeof(double) * svd_dims);

	std::cout << "Total coverage goals: " << total << std::endl;

	while (true)
	{
		int nextBest = -1;
		double bestContrib = 0;

		for (size_t i = 0; i < mtzCount(); i++)
		{
			std::vector<MtzFFTPtr>::iterator it;
			it = std::find(reordered.begin(), reordered.end(), _mtzs[i]);
			if (it != reordered.end())
			{
				continue;
			}

			double portion = 0;
			for (size_t j = 0; j < mtzCount(); j++)
			{
				wip[j] = std::max(covered[j], _clusterPtrs[i][j]);
				portion += wip[j];
			}

			double newpct = (portion / total) * 100;
			double thispct = newpct - percentage;

			if (thispct > bestContrib)
			{
				bestContrib = thispct;
				nextBest = i;
			}
		}
		
		if (nextBest == -1)
		{
			break;
		}

		double portion = 0;
		for (size_t i = 0; i < mtzCount(); i++)
		{
			covered[i] = std::max(covered[i], _clusterPtrs[nextBest][i]);
			portion += covered[i];
		}
		
		double newpct = (portion / total) * 100;
		double thispct = newpct - percentage;
		percentage = newpct;
		std::cout << "Next best contribution: " << 
		getMtzFile(nextBest)->metadata() << " with extra " << 
		thispct << "% (cumulative: " <<
		percentage << "%)." <<
		std::endl;

		reordered.push_back(_mtzs[nextBest]);

		if (_mtzs.size() == reordered.size())
		{
			break;
		}
	}

	std::cout << "Goals: " << std::endl;
	for (size_t i = 0; i < mtzCount(); i++)
	{
		std::cout << goals[i] << ", ";
	}
	std::cout << std::endl;

	std::cout << "Covered: " << std::endl;
	for (size_t i = 0; i < mtzCount(); i++)
	{
		std::cout << covered[i] << ", ";
	}
	std::cout << std::endl;

	free(wip);
	free(covered);
	
	return reordered;
}

void Group::writeHKL(std::string filename)
{
	_mySet.recip->writeHKL(filename);
}
