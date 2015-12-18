/*=====================================================================================*/
/*! \file		DistanceField.cpp
	\author		skourasm
	\brief		Implementation of class DistanceField
 */
/*=====================================================================================*/

#include "DistanceField.h"

#include <cfgMaterialUtilities.h>
using namespace cfgMaterialUtilities;

// levelset
#include "MacGrid.h"
#include "FaceField.h"
#include "LevelSet.h"
#include "Mesh.h"
#include "FAST_MARCHING_METHOD_UNIFORM.h"

const int d=3;
Type_Define_VectorD_Outside_Class(d);Type_Define_VectorDi_Outside_Class(d);Type_Define_MatrixD_Outside_Class(d);

DistanceField::DistanceField(int iDim)
{
  m_dim = iDim;
}

DistanceField::~DistanceField()
{
}

std::vector<cfgScalar> DistanceField::computeDistances(const std::vector<cfgScalar> &iPoints, std::vector<cfgScalar> *ioDerivatives)
{
  std::vector<cfgScalar> distances;

  std::vector<cfgScalar> pointsInit = iPoints;
  std::vector<cfgScalar> lengths(m_dim, 1);
  int npoint=(int)pointsInit.size()/m_dim;
  
  std::vector<cfgScalar> scalingFactors;
  bool resOk = rescaleData(pointsInit, m_dim, lengths, &scalingFactors);
  if (resOk)
  {
    Array<VectorD> points;
    for (int ipoint=0; ipoint<npoint; ipoint++)
    {
      VectorD pos;
      for (int icoord=0; icoord<d; icoord++)
      {
        pos[icoord] = pointsInit[m_dim*ipoint+icoord];
      }
      points.push_back(pos);
    }

    LevelSet<d> levelset;
    PointsToLevelSet<d> points_to_levelset(points,levelset);
    points_to_levelset.scale=32;

    points_to_levelset.Update();

    FAST_MARCHING_METHOD_UNIFORM<d> fast_marching(levelset);
    fast_marching.Fast_Marching_Method(levelset.phi);


    cfgScalar offset = (cfgScalar)levelset.grid.dx;
    //std::cout<<"dx: "<<levelset.grid.dx<<std::endl;
    //std::cout<<levelset.grid.domain_min[0]<<", "<<levelset.grid.domain_min[1]<<", "<<levelset.grid.domain_min[2]<<std::endl;
    //std::cout<<levelset.grid.domain_max[0]<<", "<<levelset.grid.domain_max[1]<<", "<<levelset.grid.domain_max[2]<<std::endl;
    for (int ipoint=0; ipoint<npoint; ipoint++)
    {
      cfgScalar distance = (cfgScalar)levelset.Phi(points[ipoint]);
      //std::cout << distance << " " << points[ipoint][0] << " " << points[ipoint][1] << " " << points[ipoint][2] << " " << std::endl;
      //distance += offset;
      distances.push_back(distance);
    }
    if (ioDerivatives)
    {
      ioDerivatives->clear();
      for (int ipoint=0; ipoint<npoint; ipoint++)
      {
        cfgScalar distance = (cfgScalar)levelset.Phi(points[ipoint]);
        VectorD normal = levelset.Normal(points[ipoint]);
        for (int icoord=0; icoord<d; icoord++)
        {
          cfgScalar val = (cfgScalar)normal[icoord];
          val *= scalingFactors[icoord];
          ioDerivatives->push_back(val);
        }
        for (int icoord=d; icoord<m_dim; icoord++)
        {
          ioDerivatives->push_back(0);
        }
      }
    }
  }
  else
  {
    distances.resize(npoint, 1);
    if (ioDerivatives)
    {
      ioDerivatives->resize(npoint*d, 0);
    }
  }
  return distances;
}



