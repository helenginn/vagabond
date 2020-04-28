#ifndef __fuck_cov__ClusterList__
#define __fuck_cov__ClusterList__

#include <vector>
#include <string>
#include <QObject>
#include <libsrc/DiffractionMTZ.h>
#include "MtzFFTPtr.h"
#include "SQLInput.h"

class QTreeWidget;
class MtzFile;
class Averager;
class QKeyEvent;
class Screen;

class ClusterList : public QObject
{
Q_OBJECT
public:
	ClusterList(QTreeWidget *widget);

	~ClusterList();

	void setCommands(std::vector<std::string> commands);
	bool loadFiles();
	void load(std::vector<DatasetPaths> paths);

	void average(Averager *item);
	void cluster(Averager *item);
	void getFromDatabase();
	void setScreen(Screen *scr)
	{
		_screen = scr;
	}
	
	void makeGroup(std::vector<MtzFFTPtr> mtzs, bool useAve);
	
	void setFiles(std::vector<std::string> files)
	{
		_filenames = files;
	}
	
	Averager *topCluster();
	void removeCluster(Averager *ave);
	void clearSelection();
	void invertSelection();
signals:
	void average();
	void cluster();
	void updateSelections();
public slots:
	void toggleDead();
	void updateColours();
	void exportAll();
	void handleResults();
	void handleError();
	void displayResults();
	void resetAverage();
	void topAverage();
	void pdbAverage();
protected:
	virtual void keyPressEvent(QKeyEvent *event);
	virtual void keyReleaseEvent(QKeyEvent *event);
private:

	Screen *_screen;
	QTreeWidget *_widget;
	QThread *_worker;
	Averager *_lastAverage;
	bool _selectMode;
	bool _removeMode;

	std::vector<std::string> _filenames;
	std::vector<Averager *> _clusters;
	std::vector<MtzFile *> _files;
	std::vector<std::string> _commands;
	double _res;
	bool _sqlInput;
};


#endif
