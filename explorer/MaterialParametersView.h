#ifndef MaterialParametersView_h
#define MaterialParametersView_h

#include <QVTKWidget.h> 
#include <vtkSmartPointer.h>
#include <QPointer>
#include <map>

#include "exProject.h"

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
class cfgPlotPoints3D;

class Resampler;

class MaterialParametersView: public QVTKWidget
{
  Q_OBJECT
public:
  MaterialParametersView();
  virtual ~MaterialParametersView();

  virtual void mouseReleaseEvent(QMouseEvent *);

  void setProject(QPointer<exProject> iProject);
  int getStructureIndex(int iPointIndex);

  void changePointsColor(const std::vector<int> &iPointIndices, int iIndPlot, const vtkVector3i &iColor);
  void changePointsColor(const std::vector<int> &iPointIndices, int iIndPlot, const std::vector<vtkVector3i> &iColors);

private:
  void highlighPoints(const std::vector<int> &iPointIndices);

  int getParameterIndex(int iStructureIndex, int iLevel);
  int getParameterDim();

  void updateChart();
  void updatePlots();

  void computeDistancesToBoundary(const std::vector<std::vector<cfgScalar> > &iParameters, int iParamDim);

  void rescaleParameters(int iIndex, float iScalingFactor, std::vector<float> &ioPoints, int iDim);
  float getMaxValue(int iIndex, const std::vector<float> &iPoints, int iDim);
  void convertToLogValues(int iIndex, const std::vector<float> &iPoints, int iDim, std::vector<float> &oPoints, float iEpsilon=0);

  void getValidParameterTypes(exProject::MicrostructureType iMicrostructureType, int iMicrostructureDimension, std::vector<exProject::MaterialParameterType> &oValidParameterTypes, std::vector<int> &oDimToRescale);
  std::string parameterTypeToString(const exProject::MaterialParameterType &iType);

  void displayPoint(int iPointIndex, int iLevel, const vtkVector3i &iColor);

  private slots:
    void onLevelVisibilityModified();
    void onLevelsModified();
    void onParamsToVisualizeModified();
    void onReducedPointModified();

    void onRadiusValueChanged(double);
    void onNbPointsValueChanged(int);

private:
  vtkContextView * m_vtkView;

  vtkSmartPointer<cfgChartXYZ> m_chart;
  std::vector<vtkSmartPointer<vtkTable> > m_tables;
  std::vector<std::vector<vtkSmartPointer<vtkPlot3D> > > m_plotsPerLevel;
  vtkSmartPointer<cfgPlotPoints3D> m_auxiliaryPlot;

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

  std::vector<std::vector<float> > m_physicalParametersPerLevel;
  std::vector<std::vector<float> > m_distancesToBoundary;
  std::vector<Resampler *> m_resamplers;
  double m_samplerRadius;
  int m_samplerNbPoints;

  std::vector<int> m_savedPointIndices;
  std::vector<vtkVector3i> m_savedColors;
  cfgPlotPoints3D * m_savedPlot;

  int m_smallDotSize;
  int m_largeDotSize;
};

#endif 


