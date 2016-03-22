#ifndef MaterialParametersView_h
#define MaterialParametersView_h

#include <QVTKWidget.h> 
#include <vtkSmartPointer.h>
#include <QPointer>
#include <map>

class vtkRenderer;
class vtkRenderView;
class vtkContextView;
class cfgPlotPoints3D;
class cfgChartXYZ;
class vtkTable;
class vtkVector3i;
class cfgPlotLine3D;
class cfgPlotSurface;
class vtkPlot3D;

class exProject;

class MaterialParametersView: public QVTKWidget
{
  Q_OBJECT
public:
  MaterialParametersView();
  virtual ~MaterialParametersView();

  virtual void mouseReleaseEvent(QMouseEvent *);

  void setProject(QPointer<exProject> iProject);
  int getStructureIndex(int iPointIndex);

private:
  vtkSmartPointer<vtkTable> createTable(const std::vector<float> &iValues1, int idim1, const std::vector<std::string> &iLabels, const std::vector<int> *iPointIndices);
  vtkSmartPointer<vtkTable> createTable(const std::vector<float> &iValues1, const std::vector<int> &iValues2, int idim1, int idim2, const std::vector<std::string> &iLabels, const std::vector<int> *iPointIndices=NULL);
  vtkSmartPointer<cfgPlotPoints3D> createPointPlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, std::vector<vtkVector3i> colors, float iWidth=5);
  vtkSmartPointer<cfgPlotPoints3D> createPointPlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, vtkVector3i &iColor, float iWidth=5);
  vtkSmartPointer<cfgPlotLine3D> createLinePlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, vtkVector3i &iColor, float iWidth=1);
  vtkSmartPointer<cfgPlotSurface> createSurfacePlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, vtkVector3i &iColor, int iAlpha=255);

  void highlighPoints(const std::vector<int> &iPointIndices);

  int getParameterIndex(int iStructureIndex, int iLevel);

  void updateChart();
  void updatePlots();

  void rescaleParameters(int iIndex, float iScalingFactor, std::vector<float> &ioPoints, int iDim);
  float getMaxValue(int iIndex, const std::vector<float> &iPoints, int iDim);

  private slots:
    void onLevelVisibilityModified();
    void onLevelsModified();
    void onParamsToVisualizeModified();

private:
  vtkContextView * m_vtkView;

  vtkSmartPointer<cfgChartXYZ> m_chart;
  vtkSmartPointer<vtkTable> m_tablePoint3D;
  std::vector<std::vector<vtkSmartPointer<vtkPlot3D> > > m_plotsPerLevel;

  std::vector<int> m_lastPointIndices;
  std::vector<int> m_plots2Levels;
  std::vector<int> m_plots2PlotIndex;
  std::vector<std::vector<std::vector<int> > > m_plotPointIndices;
  
  QPointer<exProject> m_project;

  /*std::vector<std::vector<std::vector<int> > > m_baseMaterials;
  std::vector<std::map<std::vector<int>, int> > m_baseMaterials2Indices;
  std::vector<std::vector<std::vector<int> > > m_materialAssignements;

  std::vector<float> m_physicalParameters;
  std::vector<int> m_levels;*/ 

  vtkSmartPointer<vtkTable> m_tableHighLightedPoints;
  vtkSmartPointer<cfgPlotPoints3D> m_plotHighLightedPoints;
};

#endif 


