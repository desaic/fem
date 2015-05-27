
#include "cfgPlotSurface.h"

#include <vtkContext2D.h>
#include <vtkContext3D.h>
#include <vtkPen.h>

#include <vtkObjectFactory.h>

#include "vtkUnsignedCharArray.h"

//-----------------------------------------------------------------------------
vtkStandardNewMacro(cfgPlotSurface)

//-----------------------------------------------------------------------------
cfgPlotSurface::cfgPlotSurface()
{
  this->NumberOfComponents = 0;
}

//-----------------------------------------------------------------------------
cfgPlotSurface::~cfgPlotSurface()
{
}

//-----------------------------------------------------------------------------
bool cfgPlotSurface::Paint(vtkContext2D *painter)
{
  if (!this->Visible)
  {
    return false;
  }
  // Get the 3D context.
  vtkContext3D *context = painter->GetContext3D();
  if (!context)
  {
    return false;
  }
  this->Update();

  // draw the surface
  if (this->Points.size() > 0)
  {
    context->ApplyPen(this->Pen.GetPointer());
    int npoint = (int)Points.size();

    if (this->NumberOfComponents==0)
    {
      unsigned char colRGBa[4];
      this->Pen->GetColor(colRGBa);
      colRGBa[3] = this->Pen->GetOpacity();

      this->NumberOfComponents = 4;
      this->Colors->Reset();
      this->Colors->Allocate( this->NumberOfComponents*npoint);

      int ipoint=0;
      for (ipoint=0; ipoint<npoint; ipoint++)
      {
        this->Colors->InsertNextTupleValue(&colRGBa[0]);
        this->Colors->InsertNextTupleValue(&colRGBa[1]);
        this->Colors->InsertNextTupleValue(&colRGBa[2]); 
        this->Colors->InsertNextTupleValue(&colRGBa[3]); 
      }
    }
    context->DrawTriangleMesh(this->Points[0].GetData(), npoint, this->Colors->GetPointer(0), this->NumberOfComponents);
  }
  return true;
}

