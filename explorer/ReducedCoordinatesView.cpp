#include "ReducedCoordinatesView.h" 

#include <vtkContextView.h>
#include <vtkRenderer.h>
#include <cfgChartXYZ.h>
#include <vtkContextScene.h>
#include <vtkRenderWindow.h>
#include <vtkPlot3D.h>
#include "cfgPlotPoints3D.h"

#include "QVTKWidgetUtilities.h"
using namespace QVTKWidgetUtilities;

#include "cfgMaterialUtilities.h"
#include "cfgUtilities.h"
#include "ExplorerUtilities.h"

#include "exProject.h"

#include <cfgMaterialUtilities.h>

ReducedCoordinatesView::ReducedCoordinatesView()
  :QVTKWidget()
{
  init();

  m_smallDotSize = 7;
  m_largeDotSize = 15;
}

ReducedCoordinatesView::~ReducedCoordinatesView()
{
}

void ReducedCoordinatesView::setProject(QPointer<exProject> iProject)
{
  m_project = iProject;

  connect(m_project, SIGNAL(pickedStructureModified()), this, SLOT(onSelectedPointModified()));
}

void ReducedCoordinatesView::init()
{
   GetRenderWindow();

  m_vtkView = vtkContextView::New();
  m_vtkView->SetInteractor(GetInteractor());
  SetRenderWindow(m_vtkView->GetRenderWindow());

  vtkRenderer * renderer = m_vtkView->GetRenderer();
  renderer->SetBackground(1,1,1);

  int width = this->width();
  int height = this->height();
  int size = std::min(width, height);
  int x = width/2 - size/2;
  int y = height/2 - size/2;

  vtkSmartPointer<cfgChartXYZ> chart = vtkSmartPointer<cfgChartXYZ>::New();
  chart->SetGeometry(vtkRectf(x, y, size, size));
  m_vtkView->GetScene()->AddItem(chart.GetPointer()); 
  m_chart = chart;
  updateChart();

  m_vtkView->GetRenderWindow()->SetMultiSamples(10);
  m_vtkView->GetRenderWindow()->Render();
}

void ReducedCoordinatesView::updateReducedCoordinates(const std::vector<cfgScalar> &iPoints, int iDim)
{
  m_pointIndex2SelectedPointIndex.clear();

  std::vector<std::string> labels;
  labels.push_back("x");
  labels.push_back("y");
  labels.push_back("");

  std::vector<cfgScalar> points = (iDim==3? iPoints: cfgMaterialUtilities::from2DTo3D(iPoints));
  int npoint = (int)points.size()/3;
  std::cout << "npoint = " << npoint << std::endl;
  if (iDim==2)
  {
    points[2] += 0.001; // vtk doesn't display anything if there is no variation in one dimension
  }

  if (npoint>0)
  {
    std::vector<int> & indices = m_project->getMicrostructuresIndices();
    const std::vector<std::vector<cfgScalar> > & distances = m_project->getDistancesToBoundary();

    int npoint=(int)indices.size();
    for (int ipoint=0; ipoint<npoint; ipoint++)
    {
      int indMicrostructure = indices[ipoint];
      m_pointIndex2SelectedPointIndex[indMicrostructure] = ipoint;
    }

    int indLevel = m_project->getSelectedGamutIndex();
    std::vector<vtkVector3i> allColors;
    ExplorerUtilities::getColorsForDistanceVisualization(distances[indLevel], allColors);

    std::vector<vtkVector3i> colors = cfgUtil::getSubVector(allColors, indices);
    if (colors.size()==0)
    {
      colors.resize(npoint, vtkVector3i(0, 0, 255));
    }
    ExplorerUtilities::getColorsForReducedCoordinatesVisualization(points, 3, colors);

    vtkSmartPointer<vtkTable> table =  QVTKWidgetUtilities::createTable(points, 3, labels, NULL);
    vtkSmartPointer<cfgPlotPoints3D> plot = QVTKWidgetUtilities::createPointPlot3D(table, labels[0], labels[1], labels[2], colors, 10);
    m_table = table;
    m_plot = plot;
    updateChart();
  }

  //std::vector<int> & getMicrostructuresIndices() {return m_microstructuresIndices;}
}

void ReducedCoordinatesView::updateChart()
{
  m_chart->ClearPlots();
  if (m_plot.Get() != NULL)
  {
    m_chart->AddPlot(m_plot);
  }
  if (m_auxiliaryPlot.Get()!=NULL)
  {
    m_chart->AddPlot(m_auxiliaryPlot);
  }
  update();
}

void ReducedCoordinatesView::mouseReleaseEvent(QMouseEvent * iMouseEvent)
{  
  QVTKWidget::mouseReleaseEvent(iMouseEvent);

  int pickedPointIndex = m_chart->pickedPointIndex();
  int pickedPlotIndex = m_chart->pickedPlotIndex();
  if (pickedPointIndex>=0)
  {
    m_project->setPickedReducedPoint(pickedPointIndex);
  }
}

void ReducedCoordinatesView::displayPoint(int iPointIndex, const vtkVector3i &iColor)
{
  std::cout << "ReducedCoordinatesView::displayPoint " << iPointIndex << std::endl;
  int indPlot = 0;
  cfgPlotPoints3D * plot = (cfgPlotPoints3D*)m_plot.Get();

  std::vector<int> rowIndices;
  rowIndices.push_back(iPointIndex);

  vtkSmartPointer<vtkTable> table = QVTKWidgetUtilities::createTable(plot->getTable(), rowIndices);
  m_auxiliaryPlot = createPointPlot3D(table, plot->GetXAxisLabel(), plot->GetYAxisLabel(), plot->GetZAxisLabel(), iColor, m_largeDotSize);

  updateChart();
  update();
}

void ReducedCoordinatesView::onSelectedPointModified()
{
  std::cout << "ReducedCoordinatesView::onSelectedPointModified" << std::endl;

  int indPoint = m_project->getPickedStructureIndex();
  if (indPoint>=0)
  {
    if (m_pointIndex2SelectedPointIndex.count(indPoint)>0)
    {
      int indMicrostructure = m_pointIndex2SelectedPointIndex[indPoint];
      std::cout << "indPoint = " << indPoint << "indMicrostructure = " << indMicrostructure << std::endl;

      vtkVector3i col(255, 0, 255);
      displayPoint(indMicrostructure, col);
    }
  }
}