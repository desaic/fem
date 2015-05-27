
#ifndef cfgPlotSurface_h
#define cfgPlotSurface_h

#include "vtkPlot3D.h"

class cfgPlotSurface : public vtkPlot3D
{
public:
  vtkTypeMacro(cfgPlotSurface, vtkPlot3D);
  static cfgPlotSurface * New();

  // Description:
  // Paint event for the XY plot, called whenever the chart needs to be drawn
  virtual bool Paint(vtkContext2D *painter);

protected:
  // Description:
  // This array indicates how the surface should be colored.
  vtkNew<vtkUnsignedCharArray> Colors;

protected:
  cfgPlotSurface();
  virtual ~cfgPlotSurface();
};

#endif //cfgPlotSurface_h
