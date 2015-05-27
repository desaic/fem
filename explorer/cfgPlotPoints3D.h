#ifndef cfgPlotPoints3D_h
#define cfgPlotPoints3D_h

#include <vtkPlotPoints3D.h>

class cfgPlotPoints3D: public vtkPlotPoints3D
{
public:
  vtkTypeMacro(cfgPlotPoints3D, vtkPlotPoints3D);
  static cfgPlotPoints3D * New();

  void setColor(int ipoint, vtkVector3i iColor);
  vtkVector3i getColor(int ipoint);

  void SetColors(std::vector<vtkVector3i> &iColors);

private:
 
protected:
  cfgPlotPoints3D();
  virtual ~cfgPlotPoints3D();

private:
};

#endif 






