#include "MainWindow.h"

#include <QVTKWidget.h> 
#include <vtkRenderView.h>
#include <vtkContextView.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>

#include <vtkContextScene.h>

#include <vtkChartXYZ.h>
#include <vtkFloatArray.h>
#include <vtkTable.h>
#include <vtkPlotSurface.h>
#include <vtkPlotPoints3D.h>
#include <vtkChart.h>

#include <vtkPointPicker.h>
#include <vtkPolyData.h>
#include <vtkIdTypeArray.h>

#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;

#include "MaterialParametersView.h"

#include "exProject.h"

MainWindow::MainWindow()
{
  setupUi(this);

  m_matParametersView = new MaterialParametersView();
  setCentralWidget(m_matParametersView);

  int dim = 2;
  m_project = new exProject(dim);
  //m_project->setFileDirectory("..//..//Output//");
  //m_project->setFileDirectory("..//..//Output//test_1_0//");
  //m_project->setFileDirectory("..//..//Output//test_1_0.5//");
  //m_project->setFileDirectory("..//..//Output//test_1_1//");
  // m_project->setFileDirectory("..//..//Output//test_1_0.1//");
  //m_project->setFileDirectory("..//..//Output//test_pairs_boundary//");
  //m_project->setFileDirectory("..//..//Output//test_boundaryMismatch//");

   m_project->setFileDirectory("..//..//Output//2D/");

  m_matParametersView->setProject(m_project);
  m_materialStructureView->setProject(m_project);
}

MainWindow::~MainWindow()
{
  delete m_project;
}

