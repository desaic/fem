/*=====================================================================================*/
/*! \file		DistanceField.cpp
	\author		skourasm
	\brief		Implementation of class DistanceField
 */
/*=====================================================================================*/

#include "DistanceField.h"

#include "DistanceFieldImpl.h"

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
  if (m_dim==4)
  {
    DistanceFieldImpl<3> distField(m_dim);
    distances = distField.computeDistances(iPoints, ioDerivatives);
  }
  else if (m_dim==5)
  {
    DistanceFieldImpl<4> distField(m_dim);
    distances = distField.computeDistances(iPoints, ioDerivatives);
  }
  else 
  {
    // not implemented;
    assert(0);
  }
  return distances;
}


