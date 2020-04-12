#ifndef __fuck_cov__MtzFFT__
#define __fuck_cov__MtzFFT__

#include <libsrc/DiffractionMTZ.h>
#include <libsrc/FFT.h>
#include <QTreeWidgetItem>

class MtzFile;

class MtzFFT : public VagFFT, public QTreeWidgetItem
{
public:
	MtzFFT(QTreeWidgetItem *parent, VagFFT &vag);
	MtzFFT(QTreeWidgetItem *parent, MtzFFT &vag);

	void setMtzFile(MtzFile *file)
	{
		_file = file;
	}

	MtzFile *getMtzFile()
	{
		return _file;
	}
private:
	MtzFile *_file;

};

#endif
