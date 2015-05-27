#ifndef cfgPlotLine3D_h
#define cfgPlotLine3D_h

#include "vtkPlotLine3D.h"

class cfgPlotLine3D : public vtkPlotLine3D
{
public:
  vtkTypeMacro(cfgPlotLine3D, vtkPlotLine3D);

  // Description:
  // Creates a 3D Chart object.
  static cfgPlotLine3D *New();

  // Description:
  // Paint event for the XYZ plot, called whenever the chart needs to be drawn.
  virtual bool Paint(vtkContext2D *painter);


protected:
  cfgPlotLine3D();
  ~cfgPlotLine3D();
};

#endif //cfgPlotLine3D_h
