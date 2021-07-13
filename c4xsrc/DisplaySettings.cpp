// cluster4x
// Copyright (C) 2019 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#include "DisplaySettings.h"
#include "Plot3D.h"
#include "ClusterList.h"
#include <hcsrc/FileReader.h>
#include <QCheckBox>
#include <QVBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>

DisplaySettings::DisplaySettings(QWidget *parent) : QMainWindow(parent)
{
	QVBoxLayout *vbox = new QVBoxLayout();
	QWidget *window = new QWidget(NULL);
	window->setLayout(vbox);

	{
		QCheckBox *l = new QCheckBox("Draw text:", window);
		l->setObjectName("draw_text");
		l->setChecked(Plot3D::showsText());
		vbox->addWidget(l);
	}

	{
		QCheckBox *l = new QCheckBox("Depth cue:", window);
		l->setObjectName("depth_cue");
		l->setChecked(Plot3D::usesDepth());
		vbox->addWidget(l);
	}

	{
		QHBoxLayout *hbox = new QHBoxLayout();
		QLabel *l = new QLabel("Point size:", window);
		hbox->addWidget(l);
		QLineEdit *e = new QLineEdit(window);
		float size = Plot3D::pointSize();
		e->setText(QString::fromStdString(f_to_str(size, 0)));
		e->setObjectName("point_size");
		hbox->addWidget(e);
		vbox->addLayout(hbox);
	}

	{
		QHBoxLayout *hbox = new QHBoxLayout();
		QLabel *l = new QLabel("Font size:", window);
		hbox->addWidget(l);
		QLineEdit *e = new QLineEdit(window);
		e->setPlaceholderText("8");
		float size = Plot3D::fontSize();
		e->setText(QString::fromStdString(f_to_str(size, 0)));
		e->setObjectName("font_size");
		hbox->addWidget(e);
		vbox->addLayout(hbox);
	}

	{
		QHBoxLayout *hbox = new QHBoxLayout();
		QPushButton *b = new QPushButton("Set", window);
		hbox->addWidget(b);
		vbox->addLayout(hbox);
		
		connect(b, &QPushButton::clicked, this, &DisplaySettings::run);
	}
	
	setCentralWidget(window);
}

void DisplaySettings::run()
{
	QLineEdit *e = findChild<QLineEdit *>("point_size");
	std::string text = e->text().toStdString();
	float size = atof(text.c_str());
	Plot3D::setPointSize(size);

	e = findChild<QLineEdit *>("font_size");
	text = e->text().toStdString();
	size = atof(text.c_str());
	Plot3D::setFontSize(size);

	QCheckBox *c = findChild<QCheckBox *>("depth_cue");
	bool cue = c->isChecked();
	Plot3D::setDepthCue(cue);

	c = findChild<QCheckBox *>("draw_text");
	bool draw = c->isChecked();
	Plot3D::setShowText(draw);
	
	_list->updateSelections();

	hide();
	deleteLater();
}
