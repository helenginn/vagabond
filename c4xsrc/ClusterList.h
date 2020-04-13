#ifndef __fuck_cov__ClusterList__
#define __fuck_cov__ClusterList__

#include <vector>
#include <string>
#include <QObject>
#include <libsrc/DiffractionMTZ.h>
#include "MtzFFTPtr.h"

class QTreeWidget;
class MtzFile;
class Averager;
class Screen;

class ClusterList : public QObject
{
Q_OBJECT
public:
	ClusterList(QTreeWidget *widget);

	void setCommands(std::vector<std::string> commands);
	void loadFiles(std::vector<std::string> files);

	void average(Averager *item);
	void cluster(Averager *item);
	void setScreen(Screen *scr)
	{
		_screen = scr;
	}
	
	void makeGroup(std::vector<MtzFFTPtr> mtzs, bool useAve);
	
	void removeCluster(Averager *ave);
	void clearSelection();
	void invertSelection();
signals:
	void average();
	void cluster();
	void updateSelections();
public slots:
	void exportAll();
	void handleResults();
	void handleError();
	void displayResults();
private:
	std::string findNextFilename(std::string file);

	Screen *_screen;
	QTreeWidget *_widget;
	QThread *_worker;
	Averager *_lastAverage;

	std::vector<Averager *> _clusters;
	std::vector<MtzFile *> _files;
	std::vector<std::string> _commands;
};


#endif
