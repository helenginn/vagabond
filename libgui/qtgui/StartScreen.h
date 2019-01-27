#ifndef __Vagabond__MonomerExplorer__
#define __Vagabond__MonomerExplorer__

#include <QtWidgets/qpushbutton.h>
#include <QtWidgets/qlabel.h>
#include <QtWidgets/qmainwindow.h>
#include "../../libsrc/shared_ptrs.h"
#include <QtWidgets/qlineedit.h>
#include <QtWidgets/qfiledialog.h>
#include <QtWidgets/qcheckbox.h>

class StartScreen: public QMainWindow
{
	Q_OBJECT
	
public:
	StartScreen(QWidget *parent = 0, int argc = 0, char *argv[] = NULL);
	~StartScreen();

	void finishUp();
private:
	QPushButton *_bRun;

	QPushButton *_bMtz;
	QLineEdit *_tMtz;

	QPushButton *_bModel;
	QLineEdit *_tModel;

	QPushButton *_bDir;
	QLineEdit *_tDir;
	QLineEdit *_tKick;
	QLineEdit *_tMinRes;
	QLineEdit *_tMaxRes;
	
	QCheckBox *_cRefine;
	QCheckBox *_cPosition;
	QCheckBox *_cInter;
	QCheckBox *_cIntra;
	QCheckBox *_cCbAngles;
	QCheckBox *_cCgAngles;
	QCheckBox *_cGlyAngles;
	QCheckBox *_cPeptide;
	QCheckBox *_cSidechains;
	
	OptionsPtr _options;
	QFileDialog *_fileDialogue;   
	QPushButton *_bOptionals;
	
	std::vector<QWidget *> _widgets;
	std::vector<QWidget *> _optWidgets;
	std::vector<QCheckBox *> _allToggle;
	std::vector<QCheckBox *> _angleToggle, _angle2Toggle;
	
	void makeButtons();
	bool _showOpts;
	int _argc;
	char **_argv;

	void getFile(std::string title, QString types, QLineEdit *edit);
	std::string findNewFolder();
private slots:
	void pushRun();
	void chooseMtz();
	void chooseModel();
	void chooseDir();
	void toggleOptionals();
	
	void toggleDisableOptions(int);
};

#endif
