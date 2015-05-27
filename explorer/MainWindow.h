#ifndef MainWindow_h
#define MainWindow_h

#include <QObject>
#include <QPointer>

#include "ui_explorer.h"

class MaterialParametersView;
class exProject;

class MainWindow: public QMainWindow, public Ui_MainWindow
{
  Q_OBJECT

public: 
  MainWindow();
  virtual ~MainWindow();

private:
  MaterialParametersView * m_matParametersView;

  QPointer<exProject> m_project;
};

#endif 

