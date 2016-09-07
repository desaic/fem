#ifndef ReducedCoordinatesView_h
#define ReducedCoordinatesView_h

#include <QVTKWidget.h>
#include <QPointer>

#include <vtkSmartPointer.h>
#include <cfgDefs.h>

class vtkContextView;
class cfgChartXYZ;
class vtkPlot3D;
class exProject;

class ReducedCoordinatesView: public QVTKWidget
{
  Q_OBJECT
public:
  ReducedCoordinatesView();
  virtual ~ReducedCoordinatesView();

  void setProject(QPointer<exProject> iProject);
  
  virtual void mouseReleaseEvent(QMouseEvent *);

  void updateReducedCoordinates(const std::vector<cfgScalar> &iPoints, int iDim);

private:
  void init();
  void updateChart();

private:
  vtkContextView * m_vtkView;
  vtkSmartPointer<cfgChartXYZ> m_chart;
  vtkSmartPointer<vtkPlot3D> m_plot;

  QPointer<exProject> m_project;
};

#endif 


