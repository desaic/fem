
#include "cfgPlotLine3D.h"

#include <vtkContext2D.h>
#include <vtkContext3D.h>
#include <vtkPen.h>

#include <vtkObjectFactory.h>

//-----------------------------------------------------------------------------
vtkStandardNewMacro(cfgPlotLine3D);

//-----------------------------------------------------------------------------
cfgPlotLine3D::cfgPlotLine3D()
{
}

//-----------------------------------------------------------------------------
cfgPlotLine3D::~cfgPlotLine3D()
{
}

//-----------------------------------------------------------------------------
bool cfgPlotLine3D::Paint(vtkContext2D *painter)
{
  if (!this->Visible || this->Points.size() == 0)
  {
    return false;
  }

  // Get the 3D context.
  vtkContext3D *context = painter->GetContext3D();
  if(context == NULL)
  {
    return false;
  }

  // Draw the line between the points
  context->ApplyPen(this->Pen);
  int iline=0, nline=(int)this->Points.size()/2;
  assert(nline%2==0);
  for (iline=0; iline<nline; iline++)
  {
    context->DrawPoly(this->Points[2*iline].GetData(), 2);
  }
  return this->vtkPlotPoints3D::Paint(painter);
}

