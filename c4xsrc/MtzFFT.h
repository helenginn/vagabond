#ifndef __fuck_cov__MtzFFT__
#define __fuck_cov__MtzFFT__

#include <DiffractionMTZ.h>
#include <FFT.h>
#include <QTreeWidgetItem>
#include "MtzFFTPtr.h"

class MtzFile;

class MtzFFT : public VagFFT, public QTreeWidgetItem
{
public:
	MtzFFT(QTreeWidgetItem *parent, VagFFT &vag);
	MtzFFT(QTreeWidgetItem *parent, MtzFFT &vag);
	MtzFFT(QTreeWidgetItem *parent);

	MtzFFTPtr shared_from_this()
	{
		return ToMtzFFTPtr(VagFFT::shared_from_this());
	}

	static bool isMtzFFT(QTreeWidgetItem *item)
	{
		return (dynamic_cast<MtzFFT *>(item) != NULL);
	}

	void setMtzFile(MtzFile *file)
	{
		_file = file;
	}

	MtzFile *getMtzFile()
	{
		return _file;
	}

	void updateText();
private:
	MtzFile *_file;

};

#endif
