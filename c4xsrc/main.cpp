#include <cstdlib>
#include <iostream>
#include <QApplication>
#include <QOpenGLContext>
#include <libsrc/Shouter.h>
#include "Screen.h"
#include "Group.h"
#include "commit.h"
#include "ClusterList.h"

#include <libsrc/DiffractionMTZ.h>
#include <libsrc/Options.h>

int main(int argc, char * argv[])
{
	std::cout << "Qt version: " << qVersion() << std::endl;

	QSurfaceFormat fmt;

	if (QOpenGLContext::openGLModuleType() == QOpenGLContext::LibGL) 
	{
		std::cout << "OpenGL 3.3 context" << std::endl;
		fmt.setVersion(3, 3);
		fmt.setProfile(QSurfaceFormat::CoreProfile);
	}
	else 
	{
		std::cout << "OpenGL 3.0 context" << std::endl;
		fmt.setVersion(3, 0);
	}

	std::cout << "OpenGL Version: " << fmt.version().first << "." <<
	fmt.version().second << std::endl;
	QSurfaceFormat::setDefaultFormat(fmt);

	QCoreApplication::setAttribute(Qt::AA_ShareOpenGLContexts, true);
	QApplication app(argc, argv);

	QOpenGLContext *global = QOpenGLContext::globalShareContext();
	global->setShareContext(global);
	global->create();

	setlocale(LC_NUMERIC, "C");
	srand(time(NULL));
	
	std::cout << "Vagabond version: " << VAGABOND_VERSION_COMMIT_ID << std::endl;
	
	const char hang[] = "--hang";
	OptionsPtr options = OptionsPtr(new Options(1, (const char **)&hang));
	Options::setRuntimeOptions(options);
	
	std::vector<std::string> files;
	std::vector<std::string> commands;
	
	for (int i = 1; i < argc; i++)
	{
		if (strlen(argv[i]) < 2)
		{
			continue;
		}

		if (strncmp(argv[i], "--", 2) == 0)
		{
			commands.push_back(argv[i]);
		}
		else
		{
			files.push_back(argv[i]);
		}

	}

	Screen scr(NULL);
	scr.show();
	
	ClusterList *list = scr.getList();
	list->setCommands(commands);
	list->setFiles(files);
	bool success = list->loadFiles();
	int status = 0;

	if (success)
	{
		status = app.exec();
	}
	
	return status;
}
