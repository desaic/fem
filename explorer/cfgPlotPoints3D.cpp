#include "cfgPlotPoints3D.h"

#include <vtkObjectFactory.h>
#include <vtkUnsignedCharArray.h>

vtkStandardNewMacro(cfgPlotPoints3D)

cfgPlotPoints3D::cfgPlotPoints3D()
  :vtkPlotPoints3D()
{
}

cfgPlotPoints3D::~cfgPlotPoints3D()
{
}

void cfgPlotPoints3D::setColor(int ipoint, vtkVector3i iColor)
{
  const unsigned char constRGB[3] = { iColor[0], iColor[1], iColor[2] };
  for (int ic=0; ic<3; ic++)
  {
    //this->Colors->SetTupleValue(3*ipoint+ic, &constRGB[ic]);
  }
}

vtkVector3i cfgPlotPoints3D::getColor(int ipoint)
{
  vtkVector3i color;
  for (int ic=0; ic<3; ic++)
  {
    unsigned char constRGB;
    this->Colors->GetTupleValue(3*ipoint+ic, &constRGB);
    color[ic] = (int)constRGB;
  }
  return color;
}

void cfgPlotPoints3D::SetColors(std::vector<vtkVector3i> &iColors)
{
  assert(iColors.size() == this->Points.size());

  this->NumberOfComponents = 3;
  this->Colors->Reset();

  for (unsigned int i = 0; i < this->Points.size(); ++i)
  {
    vtkVector3i & color = iColors[i];
    const unsigned char constRGB[3] = { color[0], color[1], color[2] };
    for (int ic=0; ic<3; ic++)
    {
      this->Colors->InsertNextTupleValue(&constRGB[ic]);
    }
  }
  this->Modified();
}








