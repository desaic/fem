/*=====================================================================================*/
/*! \file		ScoringFunction.cpp
	\author		skourasm
	\brief		Implementation of class ScoringFunction
 */
/*=====================================================================================*/

#include "ScoringFunction.h"
#include "DensityEstimator.h"

#include <cfgUtilities.h>
using namespace cfgUtil;

#include <cfgMaterialUtilities.h>
using namespace cfgMaterialUtilities;

ScoringFunction::ScoringFunction(int iDim)
{
  assert(iDim==3);
  m_dim = iDim;
  m_densityRadius = 0;
}

ScoringFunction::~ScoringFunction()
{
}

void ScoringFunction::setDensityRadius(cfgScalar iRadius)
{
  m_densityRadius = iRadius;
}

cfgScalar ScoringFunction::estimateDensityRadius(const std::vector<cfgScalar> &iPoints)
{
  DensityEstimator densityEstimator(m_dim);
  densityEstimator.setPoints(iPoints);
  cfgScalar radius = densityEstimator.estimateRadius();
  return radius;
}

std::vector<cfgScalar> ScoringFunction::computeDensities(const std::vector<cfgScalar> &iPoints)
{
  DensityEstimator densityEstimator(m_dim);
  densityEstimator.setPoints(iPoints);
  densityEstimator.setKernelRadius(m_densityRadius);
  densityEstimator.init();
  std::vector<cfgScalar> densities = densityEstimator.computeDensities();

  cfgScalar densityMax = *std::max_element(densities.begin(), densities.end());
  cfgScalar densityMin = *std::min_element(densities.begin(), densities.end());
  densities = mult(densities, 1/densityMax);

  return densities;
}

std::vector<cfgScalar> ScoringFunction::computePropertiesScores(const std::vector<cfgScalar> &iPoints)
{
  std::vector<cfgScalar> scores;

  int ipoint, npoint=(int)iPoints.size()/m_dim;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    Vector3S P = getVector3S(ipoint, iPoints);
    cfgScalar score = -(P[0]+P[1]+P[2]);
    scores.push_back(score);
  }
  return scores;
}

std::vector<cfgScalar> ScoringFunction::computeScores(const std::vector<cfgScalar> &iPoints)
{
  std::vector<cfgScalar> densities = computeDensities(iPoints);
  std::vector<cfgScalar> scores = computePropertiesScores(iPoints);
  int ipoint, npoint=(int)scores.size();
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    densities[ipoint] = 1-densities[ipoint];
  }

  //scores = mult(scores, densities);
  scores = densities;

  cfgScalar scoreMax = *std::max_element(scores.begin(), scores.end());
  cfgScalar scoreMin = *std::min_element(scores.begin(), scores.end());

  cfgScalar minProba = 0.01;
  
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    scores[ipoint] = (scores[ipoint]-scoreMin)/(scoreMax-scoreMin);
    scores[ipoint] *= densities[ipoint];

    scores[ipoint] = (1-minProba)*scores[ipoint] + minProba;
  }
  return scores;
}
