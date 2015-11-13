/*=====================================================================================*/
/*! \file		ScoringFunction.h
	\author		skourasm
	\brief		Declaration of class ScoringFunction
 */
/*=====================================================================================*/

#ifndef ScoringFunction_h
#define ScoringFunction_h

#include "cfgDefs.h"

class ScoringFunction 
{
public: 
  ScoringFunction(int idim);
  virtual ~ScoringFunction();

  std::vector<cfgScalar> computeScores(const std::vector<cfgScalar> &iPoints);

  void setDensityRadius(cfgScalar iRadius);
  cfgScalar estimateDensityRadius(const std::vector<cfgScalar> &iPoints);

private:
   std::vector<cfgScalar> computeDensities(const std::vector<cfgScalar> &iPoints);
   std::vector<cfgScalar> computePropertiesScores(const std::vector<cfgScalar> &iPoints);

protected:
  int m_dim;
  cfgScalar m_densityRadius;
};

#endif
