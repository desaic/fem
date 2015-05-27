#ifndef MaterialParametersView_h
#define MaterialParametersView_h

#include <QVTKWidget.h> 
#include <vtkSmartPointer.h>
#include <QPointer>

class vtkRenderer;
class vtkRenderView;
class vtkContextView;
class cfgPlotPoints3D;
class cfgChartXYZ;
class vtkTable;
class vtkVector3i;
class cfgPlotLine3D;
class cfgPlotSurface;

class exProject;

class MaterialParametersView: public QVTKWidget
{
public:
  MaterialParametersView();
  virtual ~MaterialParametersView();

  virtual void mouseReleaseEvent(QMouseEvent *);

  void setProject(QPointer<exProject> iProject);
  int getStructureIndex(int iPointIndex);

private:
  vtkSmartPointer<vtkTable> createTable(const std::vector<float> &iValues1, const std::vector<int> &iValues2, int idim1, int idim2, std::string *iLabels, const std::vector<int> *iPointIndices=NULL);
  vtkSmartPointer<cfgPlotPoints3D> createPointPlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, std::vector<vtkVector3i> colors, float iWidth=5);
  vtkSmartPointer<cfgPlotPoints3D> createPointPlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, vtkVector3i &iColor, float iWidth=5);
  vtkSmartPointer<cfgPlotLine3D> createLinePlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, vtkVector3i &iColor, float iWidth=1);
  vtkSmartPointer<cfgPlotSurface> createSurfacePlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, vtkVector3i &iColor, int iAlpha=255);

private:
  vtkContextView * m_vtkView;

  vtkSmartPointer<cfgChartXYZ> m_chart;
  vtkSmartPointer<vtkTable> m_table;
  std::vector<int> m_lastPointIndices;

  QPointer<exProject> m_project;
};

#endif 


