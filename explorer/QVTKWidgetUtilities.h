#ifndef QVTKWidgetUtilities_h
#define QVTKWidgetUtilities_h

#include <vtkSmartPointer.h>
#include <vtkVector.h>
#include <vector>

class vtkTable;
class cfgPlotPoints3D;
class cfgPlotLine3D;
class cfgPlotSurface;

namespace QVTKWidgetUtilities
{
  vtkSmartPointer<vtkTable> createTable(vtkTable *iTable, const std::vector<int> &iRowIndices);
  vtkSmartPointer<vtkTable> createTable(const std::vector<float> &iValues1, int idim1, const std::vector<std::string> &iLabels, const std::vector<int> *iPointIndices);
  vtkSmartPointer<vtkTable> createTable(const std::vector<float> &iValues1, const std::vector<int> &iValues2, int idim1, int idim2, const std::vector<std::string> &iLabels, const std::vector<int> *iPointIndices=NULL);

  vtkSmartPointer<cfgPlotPoints3D> createPointPlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, std::vector<vtkVector3i> colors, float iWidth=5);
  vtkSmartPointer<cfgPlotPoints3D> createPointPlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, const vtkVector3i &iColor, float iWidth=5);
  vtkSmartPointer<cfgPlotLine3D> createLinePlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, const vtkVector3i &iColor, float iWidth=1);
  vtkSmartPointer<cfgPlotSurface> createSurfacePlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, const vtkVector3i &iColor, int iAlpha=255);
};

#endif 


