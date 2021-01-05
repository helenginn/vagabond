#ifndef __fuck_cov__ClusterList__
#define __fuck_cov__ClusterList__

#include <vector>
#include <string>
#include <QObject>
#include <vec3.h>
#include <mat3x3.h>
#include "MtzFFTPtr.h"
#include "ExportType.h"
#include "DatasetPath.h"

class QTreeWidget;
class QTreeWidgetItem;
class AveCSV;
class MtzFile;
class Group;
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
	void loadClusters(std::string contents);
	void load(std::vector<DatasetPath> paths);
	
	void setReindexMatrix(mat3x3 reindex, vec3 translate);

	void average(Group *item);
	void cluster(Group *item);
	void getFromDatabase();
	void getFromFolders();
	void getFromStream();
	void getFromCSV(AveCSV *csv);
	void loadFromMultistate(std::string pdb);
	void setScreen(Screen *scr)
	{
		_screen = scr;
	}
	
	void makeGroup(std::vector<MtzFFTPtr> mtzs);
	
	void setFiles(std::vector<std::string> files);
	
	size_t groupCount()
	{
		return _clusters.size();
	}
	
	Group *group(int i)
	{
		return _clusters[i];
	}
	
	size_t fileCount()
	{
		return _files.size();
	}
	
	MtzFile *mtzFile(int i)
	{
		return _files[i];
	}

	Group *topCluster();
	void reorderMTZs();
	void removeCluster(Group *ave);
	void clearSelection();
	void invertSelection();
	void exportAll(ExportType type);
signals:
	void average();
	void cluster();
	void updateSelections();
public slots:
	void toggleDead();
	void switchCSV(int c);
	void cycleCSV(bool forward);
	void updateColours();
	void prepareMenu(const QPoint &p);
	void prepDirs();
	void handleResults();
	void handleError();
	void displayResults();
	void selectedResults();
	void originalAverage();
	void topAverage();
	void myAverage();

	void unitCellAverage();
	void recipAverage();
	void csvAverage();
	void pdbAverage();
protected:
	virtual void keyPressEvent(QKeyEvent *event);
	virtual void keyReleaseEvent(QKeyEvent *event);
private:
	void getFromCSV(std::string csv);

	Screen *_screen;
	QTreeWidget *_widget;
	QThread *_worker;
	Group *_lastAverage;
	bool _selectMode;
	bool _removeMode;

	std::vector<DatasetPath> _paths;
	std::vector<Group *> _clusters;
	AveCSV *_csvGroup;
	std::vector<MtzFile *> _files;
	std::vector<std::string> _commands;
	double _res;
	bool _onlyLoad;
	bool _sqlInput;
	bool _streamInput;
	bool _csvInput;
	bool _contextMenu;
	int _max;
	int _skip;
	std::string _csv;
	std::string _preload;
	std::string _stream;
	std::string _geom;
	std::string _spg;
};


#endif
