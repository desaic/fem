#include "MaterialParametersView.h"

#include <vtkRenderView.h>
#include <vtkContextView.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>

#include <vtkContextScene.h>

#include <cfgChartXYZ.h>
#include <vtkFloatArray.h>
#include <vtkTable.h>
#include <cfgPlotPoints3D.h>
#include <vtkChart.h>
#include <cfgPlotLine3D.h>
#include <cfgPlotSurface.h>

#include <vtkPointPicker.h>
#include <vtkPolyData.h>
#include <vtkIdTypeArray.h>
#include <vtkPen.h>

#include <vtkContext2D.h>
#include <vtkContext3D.h>

#include "exProject.h"

#include <cfgMaterialUtilities.h>
using namespace cfgMaterialUtilities;

#include <cfgUtilities.h>
using namespace cfgUtil;

MaterialParametersView::MaterialParametersView()
  :QVTKWidget()
{
}

MaterialParametersView::~MaterialParametersView()
{
}

vtkSmartPointer<vtkTable> MaterialParametersView::createTable(const std::vector<float> &iValues1, const std::vector<int> &iValues2, int idim1, int idim2, std::string *iLabels, const std::vector<int> *iPointIndices)
{
  vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

  assert(idim1>=3);
  int npoint = (int)(iPointIndices? iPointIndices->size(): iValues1.size()/idim1);
  int icol=0;
  for (icol=0; icol<idim1; icol++)
  {
    table->AddColumn(vtkSmartPointer<vtkFloatArray>::New());
    table->GetColumn(icol)->SetName(iLabels[icol].c_str());
  }
  for (icol=0; icol<idim2; icol++)
  {
    table->AddColumn(vtkSmartPointer<vtkIntArray>::New());
    table->GetColumn(icol+idim1)->SetName(iLabels[icol+idim1].c_str());
  }
  table->SetNumberOfRows(npoint);

  assert(iValues1.size()%idim1==0 && iValues2.size()%idim2==0 && (iValues2.size()>0?  iValues1.size()/idim1==iValues2.size()/idim2: true));
  int ipoint=0;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    int indPoint = (iPointIndices? (*iPointIndices)[ipoint]: ipoint);
    for (icol=0; icol<idim1; icol++)
    {
      table->SetValue(ipoint, icol,  iValues1[idim1*indPoint + icol]);
    }
    for (icol=0; icol<idim2;  icol++)
    {
      table->SetValue(ipoint, icol+idim1,  iValues2[idim2*indPoint + icol]);
    }
  }
  return table;
}

vtkSmartPointer<cfgPlotPoints3D> MaterialParametersView::createPointPlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, std::vector<vtkVector3i> colors, float iWidth)
{
  vtkSmartPointer<cfgPlotPoints3D> plot = vtkSmartPointer<cfgPlotPoints3D>::New();
  plot->SetInputData(iTable, iLabel1, iLabel2, iLabel3);
  plot->SetColors(colors);
  plot->GetPen()->SetWidth(iWidth);
  return plot;
}

vtkSmartPointer<cfgPlotPoints3D> MaterialParametersView::createPointPlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, vtkVector3i &iColor, float iWidth)
{
  vtkSmartPointer<cfgPlotPoints3D> plot = vtkSmartPointer<cfgPlotPoints3D>::New();
  plot->SetInputData(iTable, iLabel1, iLabel2, iLabel3);
  plot->GetPen()->SetColor(iColor[0], iColor[1], iColor[2]);
  plot->GetPen()->SetWidth(iWidth);
  return plot;
}

vtkSmartPointer<cfgPlotLine3D> MaterialParametersView::createLinePlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, vtkVector3i &iColor, float iWidth)
{
  vtkSmartPointer<cfgPlotLine3D> plot = vtkSmartPointer<cfgPlotLine3D>::New();
  plot->SetInputData(iTable, iLabel1, iLabel2, iLabel3);
  plot->GetPen()->SetColor(iColor[0], iColor[1], iColor[2]);
  plot->GetPen()->SetWidth(iWidth);
  return plot;
}

vtkSmartPointer<cfgPlotSurface> MaterialParametersView::createSurfacePlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, vtkVector3i &iColor, int iAlpha)
{
  vtkSmartPointer<cfgPlotSurface> plot = vtkSmartPointer<cfgPlotSurface>::New();
  plot->SetInputData(iTable, iLabel1, iLabel2, iLabel3);
  plot->GetPen()->SetColor(iColor[0], iColor[1], iColor[2]);
  plot->GetPen()->SetOpacity(iAlpha);
  return plot;
}

void MaterialParametersView::setProject(QPointer<exProject> iProject)
{
  m_project = iProject;

  Q_ASSERT(m_project);
  // read data
  std::vector<std::vector<Scalar> > physicalParametersPerLevel(2);
  bool ResOk = m_project->getMaterialParameters(physicalParametersPerLevel[0], 1);
  ResOk = m_project->getMaterialParameters(physicalParametersPerLevel[1], 2);

  std::vector<Scalar> physicalParameters;
  int ilevel=0, nlevel=(int)physicalParametersPerLevel.size();
  for (ilevel=0; ilevel<nlevel; ilevel++)
  {
    physicalParameters.insert(physicalParameters.end(), physicalParametersPerLevel[ilevel].begin(), physicalParametersPerLevel[ilevel].end());
  }

  std::vector<vtkVector3i> colors;
  std::vector<int> levels;

  int nparam = 3;
  std::vector<vtkVector3i> matColors;
  matColors.push_back(vtkVector3i(255, 0, 0));
  matColors.push_back(vtkVector3i(0, 0, 255));
  
  int nrow=0;
  for (ilevel=0; ilevel<nlevel; ilevel++)
  {
    nrow += (int)physicalParametersPerLevel[ilevel].size()/nparam;
    colors.resize(nrow, matColors[ilevel]);
    levels.resize(nrow, ilevel+1);
    m_lastPointIndices.push_back(nrow-1);
  }
 
  vtkSmartPointer<cfgChartXYZ> chart = vtkSmartPointer<cfgChartXYZ>::New();

  std::string labels[4]= {"Density", "Y1", "Y2", "Level"};
  vtkSmartPointer<vtkTable> table =  createTable(physicalParameters, levels, 3, 1, labels);
  vtkSmartPointer<cfgPlotPoints3D> points3D = createPointPlot3D(table, "Y1", "Y2", "Density", colors);
  chart->AddPlot(points3D);
  m_table = table;

  std::vector<int> convexHull, boundaryVertices, boundaryFaces;
  //computeBoundaryPoints(physicalParametersPerLevel[0], 3, convexHull);
  computeDelaundayTriangulation(physicalParametersPerLevel[0], 3, convexHull, boundaryVertices, boundaryFaces);

  if (boundaryVertices.size()>0)
  {
    vtkSmartPointer<vtkTable> table =  createTable(physicalParameters, levels, 3, 1, labels, &boundaryVertices);
    vtkSmartPointer<cfgPlotPoints3D> boundaryPlot = createPointPlot3D(table, "Y1", "Y2", "Density", vtkVector3i(0,255,0), 10);
    chart->AddPlot(boundaryPlot);
  }
 
  if (convexHull.size()>0)
  {
    std::vector<int> edgeIndices;
    //getEdgesFromTetFaceIndexArray(convexHull, edgeIndices);

    getEdgesFromTriFaceIndexArray(boundaryFaces, edgeIndices);

    vtkSmartPointer<vtkTable> edgesTable =  createTable(physicalParameters, levels, 3, 1, labels, &edgeIndices);
    vtkSmartPointer<cfgPlotLine3D> edgesPlot = createLinePlot3D(edgesTable, "Y1", "Y2", "Density", vtkVector3i(0,0,0));
    chart->AddPlot(edgesPlot);

    vtkSmartPointer<vtkTable> facesTable =  createTable(physicalParameters, levels, 3, 1, labels, &boundaryFaces);
    vtkSmartPointer<cfgPlotSurface> surfacePlot = createSurfacePlot3D(facesTable, "Y1", "Y2", "Density", vtkVector3i(255,0,0), 255);
    chart->AddPlot(surfacePlot);
  } 

  GetRenderWindow();
  int width = this->width();
  int height = this->height();
  int size = std::min(width, height);
  int x = width/2 - size/2;
  int y = height/2 - size/2;

  chart->SetGeometry(vtkRectf(x, y, size, size));

  m_chart = chart;
  m_vtkView = vtkContextView::New();
  m_vtkView->SetInteractor(GetInteractor());
  SetRenderWindow(m_vtkView->GetRenderWindow());

  vtkRenderer * renderer = m_vtkView->GetRenderer();
  renderer->SetBackground(1,1,1);

  m_vtkView->GetScene()->AddItem(chart.GetPointer());

  m_vtkView->GetRenderWindow()->SetMultiSamples(10);
  m_vtkView->GetRenderWindow()->Render();
}

int MaterialParametersView::getStructureIndex(int iPointIndex)
{
  int shift=0;
  int ilevel=0, nlevel=(int)m_lastPointIndices.size();
  for (ilevel=0; ilevel<nlevel; ++ilevel)
  {
    int lastPointIndex = m_lastPointIndices[ilevel];
    if (iPointIndex > lastPointIndex)
    {
      shift = m_lastPointIndices[ilevel]+1;
    }
    else
    {
      break;
    }
  }
  int indStructure = iPointIndex-shift;
  return indStructure;
}

void MaterialParametersView::mouseReleaseEvent(QMouseEvent * iMouseEvent)
{  
  QVTKWidget::mouseReleaseEvent(iMouseEvent);

  int pickedPointIndex = m_chart->pickedPointIndex();
  if (pickedPointIndex>=0)
  {
    int pickedStructureLevel = m_table->GetValue(pickedPointIndex, 3).ToInt();
    int pickedStructureIndex = getStructureIndex(pickedPointIndex);
    m_project->setPickedStructure(pickedStructureIndex, pickedStructureLevel);
  }
}



