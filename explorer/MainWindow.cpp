#include "MainWindow.h"

#include <QCheckBox.h>
#include <QDropEvent>
#include <QUrl>
#include <QDoubleSpinBox>
#include <QPushButton>

#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;

#include "MaterialParametersView.h"
#include "ReducedCoordinatesView.h"

#include "exProject.h"

MainWindow::MainWindow()
{
  setupUi(this);

  // add window for tool parameters
  QDockWidget *dockWidget_toolParams = new QDockWidget(this);
  dockWidget_toolParams->setWindowTitle("Parameters");
  QWidget *dockWidgetContent_toolParams = new QWidget();
  dockWidget_toolParams->setWidget(dockWidgetContent_toolParams);

  QGridLayout * gridLayout_toolParams = new QGridLayout(dockWidgetContent_toolParams);

  int indRow=1;
  QDoubleSpinBox * radiusSpinBox = addDoubleParameter("radius",  dockWidgetContent_toolParams, gridLayout_toolParams, indRow);
  QSpinBox * nbPointsSpinBox = addIntParameter("nb points",  dockWidgetContent_toolParams, gridLayout_toolParams, indRow);
  
  QSpacerItem * verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);
  //QPushButton * pushButton_ok = new QPushButton("OK", dockWidgetContent_toolParams);
  
  int indCol=0, rowSpan=1, colSpan=1;
  //gridLayout_toolParams->addWidget(pushButton_ok,indRow, indCol, rowSpan, 2);
  //indRow++;
  gridLayout_toolParams->addItem(verticalSpacer, indRow, indCol, rowSpan, colSpan);

  addDockWidget(static_cast<Qt::DockWidgetArea>(1), dockWidget_toolParams);
  //QMainWindow::splitDockWidget(dockWidget_toolParams, layersWidget, static_cast<Qt::Orientation>(2));

  m_matParametersView = new MaterialParametersView();
  setCentralWidget(m_matParametersView);

  m_reducedCoordinatesView = new ReducedCoordinatesView();
  QDockWidget *dockWidget = new QDockWidget(this);
  dockWidget->setWidget(m_reducedCoordinatesView);
  //dockWidget->setFeatures( QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetClosable);
  dockWidget->setObjectName(QString::fromUtf8("family plot"));
  dockWidget->setVisible(true);
  addDockWidget(static_cast<Qt::DockWidgetArea>(2), dockWidget);
  dockWidget->setFloating(true);

  setAcceptDrops(true);

  std::vector<QString> materialParameterStrings((int)exProject::UndefinedType);
  materialParameterStrings[exProject::EXType] = "EX";
  materialParameterStrings[exProject::EYType] = "EY";
  materialParameterStrings[exProject::EZType] = "EZ";
  materialParameterStrings[exProject::NuXYType] = "NuXY";
  materialParameterStrings[exProject::NuXZType] = "NuXZ";
  materialParameterStrings[exProject::NuYZType] = "NuYZ";
  materialParameterStrings[exProject::MuXYType] = "MuXY";
  materialParameterStrings[exProject::MuXZType] = "MuXZ";
  materialParameterStrings[exProject::MuYZType] = "MuYZ";
  materialParameterStrings[exProject::DensityType] = "Density";
  materialParameterStrings[exProject::StrengthType] = "Strength";

  addMaterialParameterOptions(*m_x_comboBox, materialParameterStrings);
  addMaterialParameterOptions(*m_y_comboBox, materialParameterStrings);
  addMaterialParameterOptions(*m_z_comboBox, materialParameterStrings);
  addMaterialParameterOptions(*m_col_comboBox, materialParameterStrings);

  m_action2D->setChecked(true);
  m_actionRegionSelection->setChecked(true);

  std::vector<QString> typeStrings;
  typeStrings.push_back("Cubic");
  typeStrings.push_back("Orthotropic");

  addMaterialParameterOptions(*m_type_comboBox, typeStrings);

  m_x_comboBox->setCurrentIndex((int)exProject::EXType);
  m_y_comboBox->setCurrentIndex((int)exProject::NuXYType);
  m_z_comboBox->setCurrentIndex((int)exProject::DensityType);
  m_col_comboBox->setCurrentIndex((int)exProject::MuXYType);

  m_project = new exProject();
  m_project->setFileDirectory("..//..//Output//");

  setDim();
  setType();
  setParamToVisualize(0, m_x_comboBox->currentIndex());
  setParamToVisualize(1, m_y_comboBox->currentIndex());
  setParamToVisualize(2, m_z_comboBox->currentIndex());
  setParamToVisualize(3, m_col_comboBox->currentIndex());

  m_matParametersView->setProject(m_project);
  m_materialStructureView->setProject(m_project);
  m_reducedCoordinatesView->setProject(m_project);

  connect(radiusSpinBox, SIGNAL(valueChanged(double)), m_matParametersView, SLOT(onRadiusValueChanged(double)));
  connect(nbPointsSpinBox, SIGNAL(valueChanged(int)), m_matParametersView, SLOT(onNbPointsValueChanged(int)));
}

MainWindow::~MainWindow()
{
  delete m_project;
}

QDoubleSpinBox * MainWindow::addDoubleParameter(const std::string &iLabelText,  QWidget *ioParent, QGridLayout *ioGridLayout, int &ioRowIndex)
{
  QLabel * label = new QLabel(QString::fromStdString(iLabelText), ioParent);
  QDoubleSpinBox * spinBox = new QDoubleSpinBox(ioParent);
  spinBox->setMaximum(DBL_MAX);
    
  int indCol=0, rowSpan=1, colSpan=1;
  ioGridLayout->addWidget(label, ioRowIndex, indCol, rowSpan, colSpan);
  ioGridLayout->addWidget(spinBox, ioRowIndex, indCol+1, rowSpan, colSpan);
  ioRowIndex++;

  return spinBox;
}

QSpinBox * MainWindow::addIntParameter(const std::string &iLabelText,  QWidget *ioParent, QGridLayout *ioGridLayout, int &ioRowIndex)
{
  QLabel * label = new QLabel(QString::fromStdString(iLabelText), ioParent);
  QSpinBox * spinBox = new QSpinBox(ioParent);
  spinBox->setMaximum(INT_MAX);
    
  int indCol=0, rowSpan=1, colSpan=1;
  ioGridLayout->addWidget(label, ioRowIndex, indCol, rowSpan, colSpan);
  ioGridLayout->addWidget(spinBox, ioRowIndex, indCol+1, rowSpan, colSpan);
  ioRowIndex++;

  return spinBox;
}

void MainWindow::addLevelCheckBox(QWidget * iParent, int iLevel, QString iLabel)
{
  int nboxes = (int)m_levelCheckBoxes.size();
  QCheckBox * checkBox = new QCheckBox(iParent);
  //checkBox->setChecked(Qt::Checked);
  checkBox->setChecked(Qt::Unchecked);
  checkBox->setGeometry(QRect(10, 10 + nboxes*30, 80, 20));
  if (!iLabel.isEmpty())
  {
    checkBox->setText(iLabel);
  }
  else
  {
    checkBox->setText("level" + QString::number(iLevel));
  }
  //checkBox->setProperty("idx", iLevel);
  checkBox->setProperty("idx", nboxes);
  checkBox->setChecked(true);
  checkBox->show();

  m_levelCheckBoxes.push_back(checkBox);

  connect(checkBox, SIGNAL(clicked()), this, SLOT(levelCheckBoxModified()));
}

void MainWindow::levelCheckBoxModified()
{
  QCheckBox * checkBox = dynamic_cast<QCheckBox*>(sender());
  int idx = qvariant_cast<int>(checkBox->property("idx"));
  bool isVisible = checkBox->isChecked();
  m_project->setLevelVisibility(idx, isVisible);
}

void MainWindow::dropEvent(QDropEvent *iEvent)
{
  const QMimeData* mimeData = iEvent->mimeData();
  if (mimeData->hasUrls())
  {
    QStringList pathList;
    QList<QUrl> urlList = mimeData->urls();
    int ifile=0, nfile=urlList.size();
    for (ifile=0; ifile<nfile; ifile++)
    {
      QString file = urlList.at(ifile).toLocalFile();
      QString label;
      int level=0;
      bool resOk = m_project->loadFile(file, level, label);
      if (resOk)
      {
        addLevelCheckBox(layersWidgetContents, level, label);
      }
    }
  }
}

void MainWindow::dragEnterEvent(QDragEnterEvent *iEvent)
{
  iEvent->accept();
}

void MainWindow::addMaterialParameterOptions(QComboBox &iBox, const std::vector<QString> &iStrings)
{
  int istring, nstring=(int)iStrings.size();
  for (istring=0; istring<nstring; istring++)
  {
    iBox.addItem(iStrings[istring]);
  }
}

void MainWindow::on_m_action2D_triggered()
{
  setDim();
}

void MainWindow::on_m_actionFamilyExtractor_triggered()
{
  if (m_project)
  {
    m_project->runFamilyExtractor(0);
  }
  std::vector<cfgScalar> &reducedCoords = m_project->getMicrostructuresReducedCoordinates();
  //m_reducedCoordinatesView->updateReducedCoordinates(reducedCoords, 2);
}

void MainWindow::on_m_actionFamilyVisualization_triggered()
{
  if (m_project)
  {
    m_project->runFamilyExtractor(1);
  }
  std::vector<cfgScalar> &reducedCoords = m_project->getMicrostructuresReducedCoordinates();
  m_reducedCoordinatesView->updateReducedCoordinates(reducedCoords, 2);
}

void MainWindow::on_m_type_comboBox_currentIndexChanged()
{
  setType();
}

void MainWindow::on_m_x_comboBox_currentIndexChanged()
{
  updateParamToVisualize(0, m_x_comboBox->currentIndex());
}

void MainWindow::on_m_y_comboBox_currentIndexChanged()
{
  updateParamToVisualize(1, m_y_comboBox->currentIndex());
}

void MainWindow::on_m_z_comboBox_currentIndexChanged()
{
  updateParamToVisualize(2, m_z_comboBox->currentIndex());
}

void MainWindow::on_m_col_comboBox_currentIndexChanged()
{
  updateParamToVisualize(3, m_col_comboBox->currentIndex());
}

void MainWindow::setDim()
{
  if (m_project)
  {
    int idim = (m_action2D->isChecked()? 2: 3);
    m_project->setDimension(idim);
  }
}

void MainWindow::setType()
{
  if (m_project)
  {
    int indTypeIndex = m_type_comboBox->currentIndex();
    exProject::MicrostructureType type = (indTypeIndex==0? exProject::CubicType: exProject::OrthotropicType);
    m_project->setType(type); 
  }
}

void MainWindow::setParamToVisualize(int indParam, int iParamType)
{
  if (m_project)
  {
    m_project->setParameterToVisualize(indParam, (exProject::MaterialParameterType)iParamType);
  }
}

void MainWindow::updateParamToVisualize(int indParam, int iParamType)
{
  if (m_project)
  {
    m_project->updateParameterToVisualize(indParam, (exProject::MaterialParameterType)iParamType);
  }
}

