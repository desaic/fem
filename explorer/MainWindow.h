#ifndef MainWindow_h
#define MainWindow_h

#include <QObject>
#include <QPointer>

#include "ui_explorer.h"

class QCheckBox;

class MaterialParametersView;
class exProject;

class MainWindow: public QMainWindow, public Ui_MainWindow
{
  Q_OBJECT

public: 
  MainWindow();
  virtual ~MainWindow();

protected:
  void dropEvent(QDropEvent *iEvent);
  void dragEnterEvent(QDragEnterEvent *iEvent);

private:
  void addLevelCheckBox(QWidget * iParent, int iLevel, QString iLabel="");

private slots:
  void levelCheckBoxModified();

private:
  MaterialParametersView * m_matParametersView;

  QPointer<exProject> m_project;

  std::vector<QCheckBox *> m_levelCheckBoxes;
};

#endif 

