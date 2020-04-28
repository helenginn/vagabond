#ifndef __fuck_cov__Averager__
#define __fuck_cov__Averager__

#include <libsrc/FFT.h>
#include <QTreeWidgetItem>
#include <fstream>

#include "MtzFFTPtr.h"
#include <libsrc/DiffractionMTZ.h>

typedef enum
{
	AveDiffraction,
	AveCAlpha,
} AveragerType;

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

typedef struct
{
	double sum_x;
	double sum_y;
	double sum_xx;
	double sum_yy;
	double sum_xy;
	double sum_w;
} CorrelData;

inline CorrelData empty_CD()
{
	CorrelData cd;
	memset(&cd, '\0', sizeof(cd));
	return cd;
}

inline void add_to_CD(CorrelData *cd, double &x, double &y)
{
	if (x != x || y != y)
	{
		return;
	}

	cd->sum_x += x;
	cd->sum_y += y;
	cd->sum_yy += y * y;
	cd->sum_xx += x * x;
	cd->sum_xy += x * y;
	cd->sum_w += 1;
}

inline double evaluate_CD(CorrelData &cd)
{
	double top = cd.sum_w * cd.sum_xy - cd.sum_x * cd.sum_y;
	double bottom_left = cd.sum_w * cd.sum_xx - cd.sum_x * cd.sum_x;
	double bottom_right = cd.sum_w * cd.sum_yy - cd.sum_y * cd.sum_y;
	
	double r = top / sqrt(bottom_left * bottom_right);
	
	if (r != r) return 0;
	
	return r;
}

class MatrixView;
class MtzFile;

class Averager : public QObject, public QTreeWidgetItem
{
Q_OBJECT
public:
	Averager(QTreeWidget *parent);

	void addMtz(DiffractionMtzPtr mtz, MtzFile *file, CrystalPtr c);
	void addMtz(MtzFFTPtr mtz);
	
	static bool isAverager(QTreeWidgetItem *item)
	{
		return (dynamic_cast<Averager *>(item) != NULL);
	}
	
	void setMaxResolution(double res)
	{
		_res = res;
	}
	
	void setType(AveragerType ave)
	{
		_type = ave;
	}

	AveragerType getType()
	{
		return _type;
	}

	double findPDBCorrelation(int i, int j);
	
	vec3 getCentre()
	{
		return _centre;
	}
	
	MtzFile *getMtzFile(int i);
	void setMtzSelection(size_t i, bool val);
	void flipMtzSelection(int i);
	
	MtzFFTPtr getMtz(int i)
	{
		return _mtzs[i];
	}
	
	std::vector<MtzFFTPtr> mtzs()
	{
		return _mtzs;
	}
	
	size_t mtzCount()
	{
		return _mtzs.size();
	}
	
	VagFFTPtr getAverageFFT()
	{
		return _fft;
	}
	
	void setAverageFFT(VagFFTPtr fft)
	{
		_fft = fft;
	}
	
	bool isMarked()
	{
		return _marked;
	}

	void setMarked(bool marked);
	void setDead(bool dead);
	
	void scaleIndividuals();
	
	void findIntercorrelations();
	void svd();
	void drawAxes();
	void exhaustiveCorrelations();
	void updateText();
	vec3 getPoint(int num, int a1, int a2, int a3);
	
	MatrixView *getCorrelMatrix()
	{
		return _correlMatrix;
	}
	
	double getDiagW(int i)
	{
		return _w[i];
	}
	
	void useOriginalAverage();
	void copyFromAverage(Averager *ave);
	
	double **getRawPtr()
	{
		return _origPtrs;
	}
	
	void writeToStream(std::ofstream &f, bool complete);
	
	std::string getError()
	{
		return _message;
	}

public slots:
	void performAverage();
	void performCluster();
signals:
	void resultReady();
	void failed();
private:
	void makeAverage(bool force = false);
	void makeCAlphaAverage(bool force);
	double findCorrelation(int i, int j);
	double comparePairwise(int i, int j);
	VagFFTPtr makeDifference(VagFFTPtr one, VagFFTPtr two);
	double pairwiseValue(VagFFTPtr oneTwoDiff, VagFFTPtr three, 
	                     VagFFTPtr four);
	void pairwiseValues(VagFFTPtr one, VagFFTPtr two, VagFFTPtr three, 
	                    VagFFTPtr oneTwoDiff, size_t maxVag);
	void scaleIndividualMtz(int i);
	void svdAlloc();
	void drawResults(double **data, std::string filename);
	void populatePolymer(MtzFFTPtr mtz, PolymerPtr p);

	std::vector<MtzFFTPtr> _mtzs;
	std::vector<std::string> _names;
	VagFFTPtr _fft;
	VagFFTPtr _origAve;
	
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
	
	vec3 _centre;
	
	AveragerType _type;
	AveragerType _origType;
};

#endif
