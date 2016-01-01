/*=====================================================================================*/
/*! \file		Resampler.h
	\author		skourasm
	\brief		Declaration of class Resampler
 */
/*=====================================================================================*/

#ifndef Resampler_h
#define Resampler_h

#include "cfgDefs.h"

class Resampler 
{
public: 
  Resampler();
  ~Resampler();

  void resample(const std::vector<cfgScalar> &iScores, int iNbTargetParticules, std::vector<int> &oParticules);
  void resample(cfgScalar iMinRadius, int iDim, const std::vector<cfgScalar> &iPoints, const std::vector<cfgScalar> &iScores, int iNbTargetParticules, std::vector<int> &oParticules);
  void resampleBoundary(cfgScalar iMinRadius, int iDim, const std::vector<cfgScalar> &iPoints, const std::vector<cfgScalar> &iScores, int iNbTargetParticules, std::vector<int> &oParticules);
  void resampleData(cfgScalar iMinRadius, int iDim, const std::vector<cfgScalar> &iPoints, const std::vector<cfgScalar> &iScores, std::vector<int> &oParticules);

private:
  void normalize(std::vector<cfgScalar> &ioValues);
};

#endif
