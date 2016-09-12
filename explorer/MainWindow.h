#ifndef MainWindow_h
#define MainWindow_h

#include <QObject>
#include <QPointer>

#include "ui_explorer.h"

class QCheckBox;
class QSpinBox;
class QDoubleSpinBox;

class MaterialParametersView;
class ReducedCoordinatesView;
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
  void addMaterialParameterOptions(QComboBox &iBox, const std::vector<QString> &iStrings);

  QSpinBox * addIntParameter(const std::string &iLabelText,  QWidget *ioParent, QGridLayout *ioGridLayout, int &ioRowIndex);
  QDoubleSpinBox * addDoubleParameter(const std::string &iLabelText,  QWidget *ioParent, QGridLayout *ioGridLayout, int &ioRowIndex);

  void setType();
  void setDim();
  void setParamToVisualize(int indParam, int ParamType);
  void updateParamToVisualize(int indParam, int ParamType);

private slots:
  void levelCheckBoxModified();
  void on_m_action2D_triggered();
  void on_m_actionFamilyExtractor_triggered();
  void on_m_actionFamilyVisualization_triggered();

  void on_m_type_comboBox_currentIndexChanged();
  void on_m_x_comboBox_currentIndexChanged();
  void on_m_y_comboBox_currentIndexChanged();
  void on_m_z_comboBox_currentIndexChanged();
  void on_m_col_comboBox_currentIndexChanged();

private:
  MaterialParametersView * m_matParametersView;
  ReducedCoordinatesView * m_reducedCoordinatesView;

  QPointer<exProject> m_project;

  std::vector<QCheckBox *> m_levelCheckBoxes;
};

#endif 

