#include "MainWindow.h"

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);
	app.setApplicationName("explorer");

  MainWindow mainWindow;  
	mainWindow.show();

	return app.exec();
}

