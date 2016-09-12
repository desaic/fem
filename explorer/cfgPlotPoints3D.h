#ifndef cfgPlotPoints3D_h
#define cfgPlotPoints3D_h

#include <vtkPlotPoints3D.h>

class cfgPlotPoints3D: public vtkPlotPoints3D
{
public:
  vtkTypeMacro(cfgPlotPoints3D, vtkPlotPoints3D);
  static cfgPlotPoints3D * New();

  // from vtkPlotPoints3D;
  virtual void SetInputData(vtkTable *input, const vtkStdString &xName,
                            const vtkStdString &yName,
                            const vtkStdString &zName);

  void setColor(int ipoint, vtkVector3i iColor);
  vtkVector3i getColor(int ipoint);

  void SetColors(std::vector<vtkVector3i> &iColors);

  vtkTable * getTable() {return m_table;}

private:
 
protected:
  cfgPlotPoints3D();
  virtual ~cfgPlotPoints3D();

  vtkTable * m_table;

private:
};

#endif 






