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

 /* int ilevel, nlevel=65 ;
  for (ilevel=1; ilevel<nlevel; ilevel*=2)
  {
      addLevelCheckBox(dockWidgetContents, ilevel);
  }*/ 

  int dim = 2;
  m_project = new exProject(dim);
  m_project->setFileDirectory("..//..//Output//");
  //m_project->setFileDirectory("..//..//Output//Stress_Deformation_8x8_1000_samples_Ratio10//");
  //m_project->setFileDirectory("..//..//Output//Stress_Deformation_8x8_1000_samples_Ratio1000//");

  //m_project->setFileDirectory("..//..//Output//Test1_Y1_Y2_Density_RatioStiffSoft10//");
  //m_project->setFileDirectory("..//..//Output//Test2_Nu1_Nu2_Density_RatioStiffSoft10//");
  //m_project->setFileDirectory("..//..//Output//Test3_Y1_Y2_Density_RatioStiffSoft1000//");
  //m_project->setFileDirectory("..//..//Output//Test4_Nu1_Nu2_Density_RatioStiffSoft1000//");

  //m_project->setFileDirectory("..//..//Output//test_1_0//");
  //m_project->setFileDirectory("..//..//Output//test_1_0.5//");
  //m_project->setFileDirectory("..//..//Output//test_1_1//");
  // m_project->setFileDirectory("..//..//Output//test_1_0.1//");
  //m_project->setFileDirectory("..//..//Output//test_pairs_boundary//");
  //m_project->setFileDirectory("..//..//Output//test_boundaryMismatch//");

  //m_project->setFileDirectory("..//..//Output//2D//");
  //m_project->setFileDirectory("..//..//Output//2D_Periodic//");
  //m_project->setFileDirectory("..//..//Output//2D_Periodic//1m_2m//");
  //m_project->setFileDirectory("..//..//Output//2D_Periodic//0-99999//");
  //m_project->setFileDirectory("..//..//Output//2D_Periodic//200000-299999//"); 

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
        addLevelCheckBox(dockWidgetContents, level, label);
      }
    }
  }
}

void MainWindow::dragEnterEvent(QDragEnterEvent *iEvent)
{
  iEvent->accept();
}


