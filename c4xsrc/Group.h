// cluster4x
// Copyright (C) 2019 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#ifndef __fuck_cov__Group
#define __fuck_cov__Group

#include <libsrc/FFT.h>
#include <QTreeWidgetItem>
#include <fstream>

#include "MtzFFTPtr.h"
#include <libsrc/DiffractionMTZ.h>
#include <libsrc/maths.h>

typedef enum
{
	AveDiff,
	AveCA,
	AveUC,
} GroupType;

typedef enum
{
	GroupTop,
	GroupOriginal,
	GroupMe,
} WhichGroup;

typedef struct
{
	VagFFT *one;
	VagFFT *two;
} Pair;

inline Pair makePair(VagFFTPtr one, VagFFTPtr two)
{
	Pair pair;
	pair.one = &*one; pair.two = &*two;
	return pair;
}

typedef std::map<Pair, double> PairScore;
typedef std::map<VagFFTPtr, double> VagScore;

inline bool operator<(const Pair &p, const Pair &q)
{
	if (&*p.one > &*q.one) return true;
	if (&*p.one < &*q.one) return false;
	return &*p.two > &*q.two;
	
	return true;
}

class MatrixView;
class MtzFile;
class AveDiffraction;
class AveCAlpha;
class AveUnitCell;

typedef struct
{
	AveDiffraction *recip;
	AveCAlpha *ca;
	AveUnitCell *unitCell;
} AverageSet;


class Group : public QObject, public QTreeWidgetItem
{
Q_OBJECT
public:
	Group(QTreeWidget *parent);
	
	~Group();

	void addMtz(DiffractionMtzPtr mtz, MtzFile *file, CrystalPtr c);
	void addMtz(MtzFFTPtr mtz);
	
	static bool isGroup(QTreeWidgetItem *item)
	{
		return (dynamic_cast<Group *>(item) != NULL);
	}
	
	void setMaxResolution(double res)
	{
		_res = res;
	}

	GroupType getType()
	{
		return _type;
	}

	WhichGroup getWhichGroup()
	{
		return _which;
	}

	MtzFile *getMtzFile(int i);
	void setMtzSelection(size_t i, bool val);
	void flipMtzSelection(int i);
	
	std::vector<double> getUnitCell();
	
	MtzFFTPtr getMtz(int i)
	{
		return _mtzs[i];
	}
	
	std::vector<MtzFFTPtr> mtzs()
	{
		return _mtzs;
	}
	
	void setExported(bool exp)
	{
		_exported = true;
	}
	
	bool isExported()
	{
		return _exported;
	}
	size_t mtzCount()
	{
		return _mtzs.size();
	}
	
	VagFFTPtr getAverageFFT();
	
	bool isMarked()
	{
		return _marked;
	}
	
	bool isDead()
	{
		return _dead;
	}
	
	bool isTopGroup()
	{
		return (this == _topGroup);
	}

	void setMarked(bool marked);
	void setDead(bool dead);
	
	void findIntercorrelations();
	void svd();
	void drawAxes();
	void updateText();
	vec3 getPoint(int num, int a1, int a2, int a3);
	std::vector<MtzFFTPtr> getMtzsFromSelection();
	
	MatrixView *getCorrelMatrix()
	{
		return _correlMatrix;
	}
	
	double getDiagW(int i)
	{
		return _w[i];
	}
	
	void copyFromOriginal(Group *ave);
	
	double **getRawPtr()
	{
		return _origPtrs;
	}
	
	void writeToStream(std::ofstream &f, bool complete);
	
	std::string getError()
	{
		return _message;
	}
	
	vec3 getCentre();
	AverageSet *getWorkingSet();
	void averageRs(double *rwork, double *rfree,
	               double *swork, double *sfree);
	
	void useAverageGroup(WhichGroup which)
	{
		_which = which;
	}
	
	void useAverageType(GroupType type)
	{
		_type = type;
	}

	void changeColour(double r, double b, double g, double a);
public slots:
	void performAverage();
	void performCluster();
signals:
	void resultReady();
	void failed();
private:
	void makeAverage(bool force = false);
	void calculateAllAverages(bool force = false);
	void svdAlloc();
	void drawResults(double **data, std::string filename);
	void populatePolymer(MtzFFTPtr mtz, PolymerPtr p);

	std::vector<MtzFFTPtr> _mtzs;
	std::vector<std::string> _names;
	
	AverageSet _mySet;
	AverageSet _origSet;


	
	std::map<Pair, VagScore> _vagVals;
	std::vector<vec3> _atomPos;
	std::vector<int> _atomNum;

	void cleanupSVD();
	void matAlloc(double **raw, double ***ptrs);
	
	MatrixView *_correlMatrix;

	/* maximum resolution in Angstroms */
	double _res;
	std::string _message;
	
	double **_svdPtrs;
	double *_svd;
	
	double **_origPtrs;
	double *_orig;

	double **_clusterPtrs;
	double *_cluster;

	double **_vPtrs;
	double *_v;
	
	double *_w;
	
	bool _marked;
	bool _dead;
	bool _exported;
	
	vec3 _centre;
	
	GroupType _type;
	WhichGroup _which;
	
	static Group *_topGroup;
};

#endif
