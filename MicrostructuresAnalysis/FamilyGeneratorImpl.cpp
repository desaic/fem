#include "FamilyGeneratorImpl.h"

#include <stddef.h>
#include <assert.h>
#include <fstream>
#include <iostream>

#include "cfgUtilities.h"
#include "cfgMaterialUtilities.h"

#include "MicrostructureParam2D.hpp"

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

  /*m_matAssignments.resize(2);
  m_matAssignments[0].resize(m_n[0]*m_n[1]*m_n[2], 0);
  m_matAssignments[1].resize(m_n[0]*m_n[1]*m_n[2], 1);

  cfgScalar center[3] = {0.5, 0.5, 0.5};
  cfgScalar params[3] = {0.3, 0.3, 1};
  cfgScalar thickness = 1.5/(cfgScalar)m_n[0];
  cfgScalar angle = M_PI/4;
  int mat = 1;
  addEllipse(params, center, thickness, angle ,mat, m_matAssignments[0]);*/ 
}

bool FamilyGeneratorImpl::step()
{
  bool resOk = true;

  if (m_stage == InitStage)
  {
    init();
    m_stage = MicrostructureGenerationStage;
  }
  else if (m_stage == MicrostructureGenerationStage)
  {   
    cfgScalar dx = 1/(cfgScalar)m_n[0];
    dx *= 2;

    cfgScalar square_lenght_range[2] = {2*dx, 1-2*dx};
    Vector2S rec1_length_range[2] = {Vector2S(2*dx, 2*dx), Vector2S( 1-2*dx,  1-2*dx)};
    Vector2S rec2_length_range[2] = {Vector2S(2*dx, 2*dx), Vector2S( 1-2*dx,  1-2*dx)};
    Vector2S rec3_length_range[2] = {Vector2S(2*dx, 2*dx), Vector2S( 1-2*dx,  1-2*dx)};

    int nsample = 10;
    int n = nsample-1;

    if (0)
    {
      cfgScalar squareLength = 0.4;
      Vector2S rec1_length(0.2, 0.1);
      Vector2S rec2_length(0.2, 0.5);
      Vector2S rec3_length(0.6, 0.2);

      std::vector<int> matAssignment;
      genMicrostructure_Type1(squareLength, rec1_length, rec2_length, rec3_length, matAssignment);
      bool modified = false;
      bool resOk = cfgMaterialUtilities::filterOutNonConnectedComponents(m_n[0], m_n[1], matAssignment, modified);
    }

    for (int i=0; i<nsample; i++)
    {
      std::cout << "i = " << i << std::endl;
      cfgScalar squareLength = ((n-i)*square_lenght_range[0] + i*square_lenght_range[1])/n;
      for (int j=0; j<nsample; j++)
      {
        std::cout << "j = " << j << std::endl; 
        Vector2S rec1_length;
        rec1_length[0] = ((n-j)*rec1_length_range[0][0] + j*rec1_length_range[0][1])/n;
        for (int k=0; k<nsample; k++)
        {
          rec1_length[1] = ((n-k)*rec1_length_range[1][0] + k*rec1_length_range[1][1])/n;
          for (int l=0; l<nsample; l++)
          {
            Vector2S rec2_length;
            rec2_length[0] = ((n-l)*rec2_length_range[0][0] + l*rec2_length_range[0][1])/n;
            for (int m=0; m<nsample; m++)
            {
              rec2_length[1] = ((n-m)*rec2_length_range[1][0] + m*rec2_length_range[1][1])/n;
              for (int p=0; p<nsample; p++)
              {
                Vector2S rec3_length;
                rec3_length[0] = ((n-p)*rec3_length_range[0][0] + p*rec3_length_range[0][1])/n;
                for (int q=0; q<nsample; q++)
                {
                  rec3_length[1] = ((n-q)*rec3_length_range[1][0] + q*rec3_length_range[1][1])/n;

                  std::vector<int> matAssignment;
                  genMicrostructure_Type1(squareLength, rec1_length, rec2_length, rec3_length, matAssignment);
                  bool modified = false;
                  bool resOk = cfgMaterialUtilities::filterOutNonConnectedComponents(m_n[0], m_n[1], matAssignment, modified);
                  if (resOk && !modified)
                  {
                    std::cout << "added" << std::endl;
                    m_matAssignments.push_back(matAssignment);
                  }
                }
              }
            }
          }
        }
      }
    }
    std::cout << "nmat = " << m_matAssignments.size() << std::endl;
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

// square in the middle with diagonal beams + crosses + lateral beams
void FamilyGeneratorImpl::genMicrostructure_Type1(const cfgScalar &iSquareLength, const Vector2S &iRec1_Lengths, const Vector2S &iRec2_Lengths, const Vector2S &iRec3_Lengths, std::vector<int> &oMicrostructure)
{
  int nEle = m_n[0] * m_n[1];
  oMicrostructure.clear();
  oMicrostructure.resize(nEle, 1);

  std::vector<int> arrSize;
  arrSize.push_back(m_n[0]);
  arrSize.push_back(m_n[1]);
  
  //center rectangle
  Rectangle2D r;
  r.center[0] = 0.5;
  r.center[1] = 0.5;
  r.length[0] = iSquareLength;
  r.length[1] = iSquareLength;

  //branches on diagonal of center
  Rectangle2D r1;
  r1.center[0] = 0.25;
  r1.center[1] = 0.25;
  r1.length[0] = iRec1_Lengths[0];
  r1.length[1] = iRec1_Lengths[1];
  r1.rot = M_PI/4;

  //soft beams on sides of center.
  Rectangle2D r2;
  r2.center[0] = 0.5;
  r2.center[1] = 0;
  r2.length[0] = iRec2_Lengths[0];
  r2.length[1] = iRec2_Lengths[1];

  //crosses
  Rectangle2D r3;
  r3.center[0] = 0;
  r3.center[1] = 0;
  r3.length[0] = iRec3_Lengths[0];
  r3.length[1] = iRec3_Lengths[1];

  //a test ellipse
  //Ellipse2D e;
  //e.center[0] = 0.5;
  //e.center[1] = 0.5;
  //e.r[0] = 0.4;
  //e.r[1] = 0.1;
  //e.rot = 3.14/3;
  //drawEllipse(e, arr, arrSize, 0.0);

  drawRectangle<int>(r, oMicrostructure, arrSize, 0.0);
  drawRectangle<int>(r1, oMicrostructure, arrSize, 0.0);
  drawRectangle<int>(r2, oMicrostructure, arrSize, 0.0);
  drawRectangle<int>(r3, oMicrostructure, arrSize, 0.0);

  make2DCubic(oMicrostructure, arrSize);
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






