/*=====================================================================================*/
/*! \file		Resampler.h
	\author		skourasm
	\brief		Declaration of class Resampler
 */
/*=====================================================================================*/

#ifndef Resampler_h
#define Resampler_h

#include "cfgDefs.h"
class DistanceTool;

class Resampler 
{
public: 
  Resampler();
  ~Resampler();

  void setPoints(const std::vector<cfgScalar> * iPoints, int iDim);
  void setScores(const std::vector<cfgScalar> * iScores);

  void init();
  void resampleBoundary(int iSeedIndex, cfgScalar iMinRadius, cfgScalar iMaxRadius, int iNbTargetParticules, std::vector<int> &oParticules);
  void resampleBoundary(int iSeedIndex, cfgScalar iMinRadius, cfgScalar iMaxRadius, cfgScalar iRatioDist, std::vector<int> &oParticules);

  static void resample(const std::vector<cfgScalar> &iScores, int iNbTargetParticules, std::vector<int> &oParticules, cfgScalar iRatio=0);
  static void resample(cfgScalar iMinRadius, int iDim, const std::vector<cfgScalar> &iPoints, const std::vector<cfgScalar> &iScores, int iNbTargetParticules, std::vector<int> &oParticules, cfgScalar iRatio=0);
  static void resampleBoundary(cfgScalar iMinRadius, int iDim, const std::vector<cfgScalar> &iPoints, const std::vector<cfgScalar> &iScores, int iNbTargetParticules, std::vector<int> &oParticules);
  static void resampleData(cfgScalar iMinRadius, int iDim, const std::vector<cfgScalar> &iPoints, const std::vector<cfgScalar> &iScores, std::vector<int> &oParticules);

  static void resampleUsingNormalDistribution(const std::vector<cfgScalar> &iDistances, int iNbTargetParticules, std::vector<int> &oParticules);

private:
  static void normalize(std::vector<cfgScalar> &ioValues);

private:
  const std::vector<cfgScalar> * m_scores;
  const std::vector<cfgScalar> * m_initPoints;
  std::vector<double> * m_points;
  int m_dim;

  DistanceTool * m_distanceTool;
};

#endif
