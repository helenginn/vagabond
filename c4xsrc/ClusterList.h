#ifndef __fuck_cov__ClusterList__
#define __fuck_cov__ClusterList__

#include <vector>
#include <string>
#include <QObject>
#include <hcsrc/vec3.h>
#include <hcsrc/mat3x3.h>
#include <iostream>
#include "MtzFFTPtr.h"
#include "ExportType.h"
#include "DatasetPath.h"

class QTreeWidget;
class QTreeWidgetItem;
class MyDictator;
class AveCSV;
class AveVectors;
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
	bool loadFiles(bool force = false);
	void loadClusters(std::string contents);
	void load(std::vector<DatasetPath> paths);
	
	void setReindexMatrix(mat3x3 reindex, vec3 translate);

	void addCSVSwitcher();

	void average(Group *item);
	void cluster(Group *item);
	void getFromDatabase();
	void getFromUser();
	void getFromStream();
	void getFromCSV(AveCSV *csv);
	void getFromCSV(std::string csv);
	void loadFromMultistate(std::string pdb);
	void loadFromVectorList(std::string filename);

	void setScreen(Screen *scr)
	{
		_screen = scr;
	}
	
	Group *makeGroup(std::vector<MtzFFTPtr> mtzs);
	
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
	void exportVectors();
	
	std::string valueForKey(std::string key)
	{
		if (_options.count(key) == 0)
		{
			return "";
		}

		return _options[key];
	}
	
	void addOption(std::string key, std::string value)
	{
		_options[key] = value;
		std::cout << "Setting " << key << " as " << value << std::endl;
	}
	
	void clearOptions()
	{
		_options.clear();
	}
	
	Group *getLastAverage()
	{
		return _lastAverage;
	}

	void setResolution(double res)
	{
		_res = res;
	}
	
	void setVectorList(AveVectors *v)
	{
		_vectorList = v;
	}
	
	void setDictator(MyDictator *d)
	{
		_dictator = d;
	}
	
	MyDictator *dictator()
	{
		return _dictator;
	}
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
	void exportCoordinates();
	void topAverage();
	void myAverage();

	void chooseColumns();
	void enableAllColumns();
	void autoCluster();
	void reduceDims();

	void unitCellAverage();
	void recipAverage();
	void csvAverage();
	void vecAverage();
	void pdbAverage();
protected:
	virtual void keyPressEvent(QKeyEvent *event);
	virtual void keyReleaseEvent(QKeyEvent *event);
private:

	Screen *_screen;
	QTreeWidget *_widget;
	QThread *_worker;
	Group *_lastAverage;
	bool _selectMode;
	bool _removeMode;
	std::map<std::string, std::string> _options;
	MyDictator *_dictator;

	std::vector<DatasetPath> _paths;
	std::vector<Group *> _clusters;
	AveCSV *_csvGroup;
	AveVectors *_vectorList;
	std::vector<MtzFile *> _files;
	std::vector<std::string> _commands;
	double _res;
	bool _onlyLoad;
	bool _sqlInput;
	bool _streamInput;
	bool _csvInput;
	bool _contextMenu;
	bool _mustClearSelection;

	int _max;
	int _skip;

	std::string _pdb;
	std::string _compType;
	std::string _preload;
	std::string _stream;
	std::string _geom;
	std::string _spg;
	std::string _targetID;
};


#endif
