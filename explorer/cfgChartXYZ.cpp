
#include "cfgChartXYZ.h"

#include <vtkObjectFactory.h>

#include <vtkPlot3D.h>
#include <vtkTransform.h>
#include <vtkContextMouseEvent.h>
#include <vtkVariant.h>
#include <vtkFloatArray.h>
#include <vtkContextScene.h>
#include <vtkCommand.h>
#include <vtkContextKeyEvent.h>

#include "cfgPlotPoints3D.h"

vtkStandardNewMacro(cfgChartXYZ)

cfgChartXYZ::cfgChartXYZ()
:vtkChartXYZ()
{
  m_pickingMode = false;
  m_pickedPoint = -1;
  m_pickedPlot = -1;
}

cfgChartXYZ::~cfgChartXYZ()
{
}

bool cfgChartXYZ::MouseButtonReleaseEvent(const vtkContextMouseEvent &mouse)
{
  if (m_pickingMode)
  {    
     if (m_pickedPoint>=0 && m_pickedPlot>=0)
     {
       cfgPlotPoints3D * plot = (cfgPlotPoints3D*)Plots[m_pickedPlot];
       plot->setColor(m_pickedPoint, m_savedColor);
     }
     m_pickedPoint = pickPoint(mouse, m_pickedPlot);
    if (m_pickedPoint>=0 && m_pickedPlot>=0)
    {
      cfgPlotPoints3D * plot = (cfgPlotPoints3D*)Plots[m_pickedPlot];
      m_savedColor = plot->getColor(m_pickedPoint);
      plot->setColor(m_pickedPoint, vtkVector3i(0, 255, 0));
    } 
    this->Scene->SetDirty(true);
    this->InvokeEvent(vtkCommand::InteractionEvent);
  }
  return vtkChartXYZ::MouseButtonReleaseEvent(mouse);
}

bool cfgChartXYZ::MouseButtonPressEvent(const vtkContextMouseEvent &mouse)
{
   m_pickingMode = false;
  if (mouse.GetButton() == vtkContextMouseEvent::LEFT_BUTTON && mouse.GetModifiers() == vtkContextMouseEvent::CONTROL_MODIFIER)
  {
     m_pickingMode = true;
    return true;
  }
  return false;
}

bool cfgChartXYZ::MouseMoveEvent(const vtkContextMouseEvent &mouse)
{
  if (mouse.GetButton() == vtkContextMouseEvent::LEFT_BUTTON)
  {
    if (mouse.GetModifiers() == vtkContextMouseEvent::CONTROL_MODIFIER)
    {
      //return this->Rotate(mouse);
    }
    else
    {
      //return this->Spin(mouse);
      return this->Rotate(mouse);
    }
  }
  if (mouse.GetButton() == vtkContextMouseEvent::MIDDLE_BUTTON)
  {
    //if (mouse.GetModifiers() == vtkContextMouseEvent::SHIFT_MODIFIER)
    {
      return this->Pan(mouse);
    }
    //else
    {
     // return this->Zoom(mouse);
    }
  }
  return false;
} 

int cfgChartXYZ::pickPoint(const vtkContextMouseEvent &mouse, int &oPickedPlot)
{
  vtkVector2f screenPos(mouse.GetScreenPos().Cast<float>().GetData());

  float minLength = VTK_FLOAT_MAX;
  int closestPointIndex = -1;
  oPickedPlot = -1;

  vtkTransform * transfo = this->ContextTransform.GetPointer();
  for (unsigned int i = 0; i < this->Plots.size(); ++i)
  {
    vtkPlot3D * plot = this->Plots[i];
    std::vector<vtkVector3f> points = plot->GetPoints();
    int ipoint=0, npoint=(int)points.size(); 
    for (ipoint=0; ipoint<npoint; ipoint++)
    {
      vtkVector3f & p = points[ipoint];
      float coord[3];
      this->ContextTransform->TransformPoint(p.GetData(), coord);
      vtkVector3f q(coord);

      vtkVector2f u(q[0]-screenPos[0], q[1]-screenPos[1]);
      float l = u.SquaredNorm();
      if (l < minLength)
      {
        minLength = l;
        closestPointIndex = ipoint;
        oPickedPlot = i;
      }
    }
  }
  return closestPointIndex;
}

bool cfgChartXYZ::KeyPressEvent(const vtkContextKeyEvent &key)
{
  if (0) //key.GetKeyCode()=='x')
  {
    LookDownX();
  }
  else
  {
    return  vtkChartXYZ::KeyPressEvent(key);
  }
  return true;
}

void cfgChartXYZ::LookDownX()
{
  this->InvokeEvent(vtkCommand::InteractionEvent);
  this->Rotation->Identity();
  this->Rotation->RotateX(-90.0);
  //this->Rotation->RotateZ(-45.0);
  this->Scene->SetDirty(true);
}
