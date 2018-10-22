//
//  QTinker.h
//  Windexing
//
//  Created by Helen Ginn on 11/12/2016.
//  Copyright (c) 2017 Helen Ginn. All rights reserved.
//

#ifndef __Windexing__QUnitCell__
#define __Windexing__QUnitCell__

#include <stdio.h>

#include <QtWidgets/qwidget.h>
#include <QtCore/qglobal.h>
#include <QtWidgets/qapplication.h>
#include <QtWidgets/qpushbutton.h>
#include <QtWidgets/qlabel.h>
#include <QtWidgets/qlineedit.h>
#include <QtWidgets/qmainwindow.h>

class VagWindow;

typedef enum
{
	DialogueUndefined,
	DialogueBMultiplier,
	DialogueRefit,
} DialogueType;

typedef void (*SimpleSet)(double);

class Dialogue : public QMainWindow
{
	Q_OBJECT

public:
	Dialogue(QWidget *parent = 0, std::string windowText = "",
	         std::string labelText = "", std::string defaultText = "",
	         std::string buttonText = "");
	QPushButton *bDialogue;
	QLineEdit *tDialogue;
	QLabel *lDialogue;

	~Dialogue();

	void setWindow(VagWindow *window)
	{
		_window = window;
	}

	void setTag(DialogueType type)
	{
		_type = type;
	}

	void setFunction(void (*func)(double))
	{
		_func = func;
	}
	
	SimpleSet getFunction()
	{
		return _func;
	}
private slots:
	void returnClicked();

private:
	VagWindow *_window;
	DialogueType _type;
	SimpleSet _func;
};

#endif /* defined(__CaroCode__QTinker__) */
