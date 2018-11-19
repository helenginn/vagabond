#ifndef __Vagabond__MonomerExplorer__
#define __Vagabond__MonomerExplorer__

#include <QtWidgets/qpushbutton.h>
#include <QtWidgets/qlabel.h>
#include <QtWidgets/qmainwindow.h>
#include "../../libsrc/shared_ptrs.h"
#include <QtWidgets/qlineedit.h>
#include <QtWidgets/qfiledialog.h>

class StartScreen: public QMainWindow
{
	Q_OBJECT
	
public:
	StartScreen(QWidget *parent = 0, int argc = 0, char *argv[] = NULL);
	~StartScreen();
private:
	QPushButton *_bRun;

	QPushButton *_bMtz;
	QLineEdit *_tMtz;
	QLabel *_lMtz;

	QPushButton *_bModel;
	QLineEdit *_tModel;
	QLabel *_lModel;

	QPushButton *_bDir;
	QLineEdit *_tDir;
	QLabel *_lDir;

	QLabel *_lTweakables;
	QLabel *_lInputs;

	QLineEdit *_tKick;
	QLabel *_lKick;
	QLabel *_lKickTip;
	
	QLineEdit *_tMinRes;
	QLabel *_lMinRes;
	
	QLineEdit *_tMaxRes;
	QLabel *_lMaxRes;
	
	OptionsPtr _options;
	QFileDialog *_fileDialogue;   
	
	void makeButtons();
	int _argc;
	char **_argv;

	void getFile(std::string title, QString types, QLineEdit *edit);
	std::string findNewFolder();
private slots:
	void pushRun();
	void chooseMtz();
	void chooseModel();
	void chooseDir();
};

#endif
