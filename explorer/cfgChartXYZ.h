#ifndef cfgChartXYZ_h
#define cfgChartXYZ_h

#include <vtkChartXYZ.h>

class cfgChartXYZ: public vtkChartXYZ
{
public:
  vtkTypeMacro(cfgChartXYZ, vtkChartXYZ);

  static cfgChartXYZ * New();

  // from vtkChartXYZ
  virtual bool MouseButtonPressEvent(const vtkContextMouseEvent &mouse);
  virtual bool MouseMoveEvent(const vtkContextMouseEvent &mouse);
  virtual bool MouseButtonReleaseEvent(const vtkContextMouseEvent &mouse);

  int pickedPointIndex() {return m_pickedPoint;}
  int pickedPlotIndex() {return m_pickedPlot;}

  
// Description:
  // Key press event.  This allows the user to snap the chart to one of three
  // different 2D views.  "x" changes the view so we're looking down the X axis.
  // Similar behavior occurs for "y" or "z".
  virtual bool KeyPressEvent(const vtkContextKeyEvent &key);

protected:
  // Description:
  // Adjust the rotation of the chart so that we are looking down the X axis.
  void LookDownX();

private:
  int pickPoint(const vtkContextMouseEvent &mouse, int &oPickedPlot);

protected:
  cfgChartXYZ();
  ~cfgChartXYZ();

private:
  bool m_pickingMode;
  int m_pickedPoint;
  int m_pickedPlot;
  vtkVector3i m_savedColor;
};

#endif 






