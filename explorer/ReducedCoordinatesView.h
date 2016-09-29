#ifndef ReducedCoordinatesView_h
#define ReducedCoordinatesView_h

#include <QVTKWidget.h>
#include <QPointer>

#include <vtkSmartPointer.h>
#include <cfgDefs.h>
#include <vtkVector.h>

class vtkContextView;
class cfgChartXYZ;
class vtkPlot3D;
class exProject;
class cfgPlotPoints3D;
class vtkTable;

class ReducedCoordinatesView: public QVTKWidget
{
  Q_OBJECT
public:
  ReducedCoordinatesView();
  virtual ~ReducedCoordinatesView();

  void setProject(QPointer<exProject> iProject);
  
  virtual void mouseReleaseEvent(QMouseEvent *);

  void updateReducedCoordinates(const std::vector<cfgScalar> &iPoints, int iDim);

  void displayPoint(int iPointIndex, const vtkVector3i &iColor);

private:
  void init();
  void updateChart();

private slots:
  void onSelectedPointModified();

private:
  vtkContextView * m_vtkView;
  vtkSmartPointer<cfgChartXYZ> m_chart;
  vtkSmartPointer<vtkPlot3D> m_plot;
  vtkSmartPointer<vtkTable> m_table;
  vtkSmartPointer<cfgPlotPoints3D> m_auxiliaryPlot;
  std::map<int, int> m_pointIndex2SelectedPointIndex;

  QPointer<exProject> m_project;

  int m_smallDotSize;
  int m_largeDotSize;
};

#endif 


