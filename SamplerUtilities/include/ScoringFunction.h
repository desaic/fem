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

  void setCacheDensities(bool iCache);
  void setNewPointIndexStart(int iNewPointsIndex);
  std::vector<cfgScalar> computeScores(const std::vector<cfgScalar> &iPoints);

  void setUseDistanceField(bool iUseField);
  void setDensityRadius(cfgScalar iRadius);
  cfgScalar estimateDensityRadius(const std::vector<cfgScalar> &iPoints);

private:
   std::vector<cfgScalar> computeDensities(const std::vector<cfgScalar> &iPoints);
   std::vector<cfgScalar> computePropertiesScores(const std::vector<cfgScalar> &iPoints);

protected:
  int m_dim;
  cfgScalar m_densityRadius;
  bool m_useDistanceField;
  bool m_cacheDensities;
  std::vector<cfgScalar> m_cachedUnscaledDensities;
  int m_newPointIndexStart;
};

#endif
