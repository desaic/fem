#include "MainWindow.h"

#include <QCheckBox.h>
#include <QDropEvent>
#include <QUrl>

#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;

#include "MaterialParametersView.h"

#include "exProject.h"

MainWindow::MainWindow()
{
  setupUi(this);

  m_matParametersView = new MaterialParametersView();
  setCentralWidget(m_matParametersView);

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
}

MainWindow::~MainWindow()
{
  delete m_project;
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

