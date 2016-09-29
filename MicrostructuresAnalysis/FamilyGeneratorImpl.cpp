#include "FamilyGeneratorImpl.h"

#include <stddef.h>
#include <assert.h>
#include <fstream>
#include <iostream>

#include "cfgUtilities.h"
#include "cfgMaterialUtilities.h"

FamilyGeneratorImpl::FamilyGeneratorImpl()
{
  m_stage = InitStage;

  m_dim = 2;

  m_n[0] = 0;
  m_n[1] = 0;
  m_n[2] = 1;
}

FamilyGeneratorImpl::~FamilyGeneratorImpl()
{
}

void FamilyGeneratorImpl::setMicrostructureSize(int iDim, int n[3])
{
  m_dim = iDim;
  for (int idim=0; idim<m_dim; idim++)
  {
    m_n[idim] = n[idim];
  }
}

void FamilyGeneratorImpl::init()
{
  m_matAssignments.clear();
  m_parameters.clear();

  m_matAssignments.resize(2);
  m_matAssignments[0].resize(m_n[0]*m_n[1]*m_n[2], 0);
  m_matAssignments[1].resize(m_n[0]*m_n[1]*m_n[2], 1);
  m_parameters.resize(4, 0.001);
  m_parameters.resize(8, 1);

  cfgScalar center[3] = {0.5, 0.5, 0.5};
  cfgScalar params[3] = {0.3, 0.3, 1};
  cfgScalar thickness = 1.5/(cfgScalar)m_n[0];
  cfgScalar angle = M_PI/4;
  int mat = 1;
  addEllipse(params, center, thickness, angle ,mat, m_matAssignments[0]);
}

bool FamilyGeneratorImpl::step()
{
  bool resOk = true;

  if (m_stage == InitStage)
  {
    init();
    m_stage = EndStage;
  }

  return resOk;
}

bool FamilyGeneratorImpl::run()
{
  bool resOk = true;
  while (m_stage != EndStage && resOk)
  {
    resOk = step();
  }
  return resOk;
}

// x^2/a^2 + y^2/b^2 + z^2/c^2= 1 params = [a b c]
void FamilyGeneratorImpl::addEllipse(cfgScalar *iParams, cfgScalar *iCenterCoords, cfgScalar iThickness, cfgScalar iRotationAngle, int iMaterial, std::vector<int> &ioMatAssignment)
{
  int ind = 0;
  int ijk[3];
  for (ijk[0]=0; ijk[0]<m_n[0]; ijk[0]++)
  {
    for (ijk[1]=0; ijk[1]<m_n[1]; ijk[1]++)
    {
      for (ijk[2]=0; ijk[2]<m_n[2]; ijk[2]++)
      {
        cfgScalar fmin=0, fmax=0;
        
        cfgScalar x[3];
        for (int idim=0; idim<m_dim; idim++)
        {
          x[idim] = (ijk[idim]+0.5)/(cfgScalar)m_n[idim] - iCenterCoords[idim];
        }
        if (iRotationAngle != 0)
        {
          cfgScalar X = cos(iRotationAngle)*x[0] - sin(iRotationAngle)*x[1];
          cfgScalar Y = sin(iRotationAngle)*x[0] + cos(iRotationAngle)*x[1];
          x[0] = X;
          x[1] = Y;
        }
        
        for (int idim=0; idim<m_dim; idim++)
        {
          cfgScalar rmin2 = iParams[idim]-iThickness/2;
          rmin2*=rmin2;
          cfgScalar rmax2 = iParams[idim]+iThickness/2;
          rmax2*=rmax2;

          if (idim==0)
          {
            fmin += x[idim]*x[idim]/rmin2;
            fmax += x[idim]*x[idim]/rmax2;
          }
          else
          {
            fmin -= x[idim]*x[idim]/rmin2;
            fmax -= x[idim]*x[idim]/rmax2;
          }
        }
        fmin -= 1;
        fmax -= 1;
        if (fmin>=0 && fmax<=0)
        {
          ioMatAssignment[ind] = iMaterial;
        }
        ind++;
      }
    }
  }
}






