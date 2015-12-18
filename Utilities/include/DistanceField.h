/*=====================================================================================*/
/*! \file		DistanceField.h
	\author		skourasm
	\brief		Declaration of class DistanceField
 */
/*=====================================================================================*/

#ifndef DistanceField_h
#define DistanceField_h

#include "cfgDefs.h"

class DistanceField 
{
public: 
  DistanceField(int idim);
  virtual ~DistanceField();

  std::vector<cfgScalar> computeDistances(const std::vector<cfgScalar> &iPoints, std::vector<cfgScalar> *ioDerivatives=NULL);

protected:
  int m_dim;
};

#endif
