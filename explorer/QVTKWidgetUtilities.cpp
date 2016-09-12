#include "QVTKWidgetUtilities.h"

#include <vtkTable.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>

#include "cfgPlotPoints3D.h"
#include "cfgPlotLine3D.h"
#include "cfgPlotSurface.h"
#include "vtkPen.h"

vtkSmartPointer<vtkTable> QVTKWidgetUtilities::createTable(const std::vector<float> &iValues1, int idim1, const std::vector<std::string> &iLabels, const std::vector<int> *iPointIndices)
{
  vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

  int npoint = (int)(iPointIndices? iPointIndices->size(): iValues1.size()/idim1);
  int icol=0;
  for (icol=0; icol<idim1; icol++)
  {
    table->AddColumn(vtkSmartPointer<vtkFloatArray>::New());
    table->GetColumn(icol)->SetName(iLabels[icol].c_str());
  }
  table->SetNumberOfRows(npoint);

  assert(iValues1.size()%idim1==0);
  int ipoint=0;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    int indPoint = (iPointIndices? (*iPointIndices)[ipoint]: ipoint);
    for (icol=0; icol<idim1; icol++)
    {
      table->SetValue(ipoint, icol,  iValues1[idim1*indPoint + icol]);
    }
  }
  return table;
}

vtkSmartPointer<vtkTable> QVTKWidgetUtilities::createTable(vtkTable *iTable, const std::vector<int> &iRowIndices)
{
  assert(iTable);

  vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
   
  int ncol = iTable->GetNumberOfColumns();
  for (int icol=0; icol<ncol; icol++)
  {
    table->AddColumn(vtkSmartPointer<vtkFloatArray>::New());
    table->GetColumn(icol)->SetName(iTable->GetColumnName(icol));
  }
  int nrow = iRowIndices.size();
  table->SetNumberOfRows(nrow);

  for (int irow=0; irow<nrow; irow++)
  {
    int indRow = iRowIndices[irow];
    for (int icol=0; icol<ncol; icol++)
    {
      table->SetValue(irow, icol,  iTable->GetValue(indRow, icol));
    }
  }
  return table;
}

vtkSmartPointer<vtkTable> QVTKWidgetUtilities::createTable(const std::vector<float> &iValues1, const std::vector<int> &iValues2, int idim1, int idim2, const std::vector<std::string> &iLabels, const std::vector<int> *iPointIndices)
{
  vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

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

vtkSmartPointer<cfgPlotPoints3D> QVTKWidgetUtilities::createPointPlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, std::vector<vtkVector3i> colors, float iWidth)
{
  vtkSmartPointer<cfgPlotPoints3D> plot = vtkSmartPointer<cfgPlotPoints3D>::New();
  plot->SetInputData(iTable, iLabel1, iLabel2, iLabel3);
  plot->SetColors(colors);
  plot->GetPen()->SetWidth(iWidth);
  return plot;
}

vtkSmartPointer<cfgPlotPoints3D> QVTKWidgetUtilities::createPointPlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, const vtkVector3i &iColor, float iWidth)
{
  vtkSmartPointer<cfgPlotPoints3D> plot = vtkSmartPointer<cfgPlotPoints3D>::New();
  plot->SetInputData(iTable, iLabel1, iLabel2, iLabel3);
  plot->GetPen()->SetColor(iColor[0], iColor[1], iColor[2]);
  plot->GetPen()->SetWidth(iWidth);
  return plot;
}

vtkSmartPointer<cfgPlotLine3D> QVTKWidgetUtilities::createLinePlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, const vtkVector3i &iColor, float iWidth)
{
  vtkSmartPointer<cfgPlotLine3D> plot = vtkSmartPointer<cfgPlotLine3D>::New();
  plot->SetInputData(iTable, iLabel1, iLabel2, iLabel3);
  plot->GetPen()->SetColor(iColor[0], iColor[1], iColor[2]);
  plot->GetPen()->SetWidth(iWidth);
  return plot;
}

vtkSmartPointer<cfgPlotSurface> QVTKWidgetUtilities::createSurfacePlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, const vtkVector3i &iColor, int iAlpha)
{
  vtkSmartPointer<cfgPlotSurface> plot = vtkSmartPointer<cfgPlotSurface>::New();
  plot->SetInputData(iTable, iLabel1, iLabel2, iLabel3);
  plot->GetPen()->SetColor(iColor[0], iColor[1], iColor[2]);
  plot->GetPen()->SetOpacity(iAlpha);
  return plot;
}

