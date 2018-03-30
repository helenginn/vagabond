#ifndef __Vagabond__MonomerExplorer__
#define __Vagabond__MonomerExplorer__

#include <QtWidgets/qlabel.h>
#include <QtWidgets/qlistwidget.h>
#include "../../libsrc/shared_ptrs.h"
#include "AtomListItem.h"
#include "SetterEdit.h"
#include <QtWidgets/qpushbutton.h>
#include <QtWidgets/qslider.h>
#include "../GLKeeper.h"
#include <map>

typedef struct
{
	ParamOptionType optionType;
	int scale;
	double value;
	int isZero;
	QLabel *lOpt;
	QLabel *lVal;
	const char *unit;
} ParamOption;

typedef std::map<QSlider *, ParamOption> OptionMap;

class MonomerExplorer : public QWidget
{
	Q_OBJECT

public:
	MonomerExplorer(QWidget *parent = NULL, MonomerPtr monomer = MonomerPtr())
	: QWidget(parent)
	{
		initialise(monomer);
	}

	MonomerPtr getMonomer()
	{
		return _monomer;
	}
	GLKeeper *getKeeper()
	{
		return _keeper;
	}

	void setKeeper(GLKeeper *keeper);
	~MonomerExplorer();
private slots:
	void clickedAtomListItem();
	void pushRefineDensity();
	void pushSidechainsToEnd();
	void pushSqueezeToEnd();
	void pushModelPosToEnd();
	void pushRefineToEnd();
	void setSliderValue();
private:
	void initialise(MonomerPtr monomer);
	void populateList();
	void makeRefinementButtons();
	void makeSlider(ParamOptionType option, int num, QString name, int min, int max, int scale, int defVal, const char *unit);
	void applyParamOptions(SamplerPtr sampled);
	Notifiable *preparePolymer();

	MonomerPtr _monomer;
	GLKeeper *_keeper;

	QListWidget *_atomList;

	QLabel *_lModel;
	QLabel *_lTorsion;
	QLabel *_lKick;
	QLabel *_lDampen;
	QLabel *_lPhi;
	QLabel *_lPsi;
	QLabel *_lRefineOpts;

	SetterEdit *_tTorsion;
	SetterEdit *_tKick;
	SetterEdit *_tDampen;
	SetterEdit *_tPhi;
	SetterEdit *_tPsi;

	QPushButton *_bRefineDensity;
	QPushButton *_bRefineToEnd;
	QPushButton *_bSidechainsToEnd;
	QPushButton *_bSqueezeToEnd;
	QPushButton *_bModelPosToEnd;
	QLabel *_lCorrel;

	OptionMap _optionMap;
};


#endif 
