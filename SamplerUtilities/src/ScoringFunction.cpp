/*=====================================================================================*/
/*! \file		ScoringFunction.cpp
	\author		skourasm
	\brief		Implementation of class ScoringFunction
 */
/*=====================================================================================*/

#include "ScoringFunction.h"
#include "DensityEstimator.h"
#include "DistanceField.h"

#include <cfgUtilities.h>
using namespace cfgUtil;

#include <cfgMaterialUtilities.h>
using namespace cfgMaterialUtilities;

ScoringFunction::ScoringFunction(int iDim)
{
  m_dim = iDim;
  m_densityRadius = 0;
  m_useDistanceField = 0;
  m_cacheDensities = false;
  m_newPointIndexStart = 0;
}

ScoringFunction::~ScoringFunction()
{
}

void ScoringFunction::setCacheDensities(bool iCache)
{
  m_cacheDensities = iCache;
}
 
void ScoringFunction::setNewPointIndexStart(int iNewPointIndex)
{
  m_newPointIndexStart = iNewPointIndex;
}

void ScoringFunction::setUseDistanceField(bool iUseField)
{
  m_useDistanceField = iUseField;
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
  std::vector<cfgScalar> densities;
  if (m_cacheDensities)
  {
    densities = m_cachedUnscaledDensities;
    std::vector<cfgScalar> newDensities = densityEstimator.computeDensities(m_newPointIndexStart, (int)iPoints.size()/m_dim, densities);
    densities.insert(densities.end(), newDensities.begin(), newDensities.end());
  }
  else
  {
    densities = densityEstimator.computeDensities();
  }
  if (m_cacheDensities)
  {
    m_cachedUnscaledDensities = densities;
  }
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
  /*  cfgScalar score = -(iPoints[m_dim*ipoint]+iPoints[m_dim*ipoint+1]+iPoints[m_dim*ipoint+2]);
    if (m_dim==4)
      score -= iPoints[m_dim*ipoint+3];
    scores.push_back(score);
    */ 
    // Favors materials with lower density but higher Young modulus
    cfgScalar score = -iPoints[m_dim*ipoint]+iPoints[m_dim*ipoint+1]+iPoints[m_dim*ipoint+2];
    if (m_dim==4)
      score += iPoints[m_dim*ipoint+3];
    scores.push_back(score);
  }
  return scores;
}

std::vector<cfgScalar> ScoringFunction::computeScores(const std::vector<cfgScalar> &iPoints)
{
  std::vector<cfgScalar> scores;
  if (m_useDistanceField)
  {
    std::vector<cfgScalar> points = iPoints;
    std::vector<cfgScalar> lengths(m_dim, 1);
    std::vector<int> dimToScale;
    if (1)
    {
      if (m_dim==4)
      {
        dimToScale.push_back(1);
      }
      else if (m_dim==5)
      {
        dimToScale.push_back(1);
        dimToScale.push_back(2);
      }
      else
      {
        dimToScale.push_back(1);
        dimToScale.push_back(2);
        dimToScale.push_back(3);
      }
    }
    cfgScalar eps = 1.e-6;
    convertToLogValues(points, m_dim, dimToScale, eps);
    rescaleData(points, m_dim, lengths);
    if (m_densityRadius==0)
    {
      m_densityRadius = estimateDensityRadius(points);
    }
    
    std::cout << "computing distances..." << std::endl;
    DistanceField distanceField(m_dim);
    scores = distanceField.computeDistances(iPoints);

    bool useDensity = true;
    std::vector<cfgScalar> densities;
    if (useDensity)
    {
      std::cout << "computing densities..." << std::endl;
      densities = computeDensities(points);
      cfgScalar scoreMax = *std::max_element(scores.begin(), scores.end());
      cfgScalar scoreMin = *std::min_element(scores.begin(), scores.end());

      cfgScalar minProba = 0 ;//0.1;

      if (scoreMax-scoreMin > 1.e-6)
      {
        int ipoint, npoint=(int)scores.size();
        for (ipoint=0; ipoint<npoint; ipoint++)
        {
          scores[ipoint] = (scores[ipoint]-scoreMin)/(scoreMax-scoreMin);
          scores[ipoint] /= (densities[ipoint]*densities[ipoint]);
          //scores[ipoint] /= densities[ipoint];
          scores[ipoint] = (1-minProba)*scores[ipoint] + minProba;
          //scores[ipoint] *= scores[ipoint];
        }
      }
    }
  }
  else
  {
    std::vector<cfgScalar> points = iPoints;
    std::vector<cfgScalar> lengths(m_dim, 1);
    rescaleData(points, m_dim, lengths);

    if (m_densityRadius==0)
    {
      m_densityRadius = estimateDensityRadius(points);
    }
    std::cout << "density radius = " << m_densityRadius << std::endl;
    std::vector<cfgScalar> densities = computeDensities(points);
    cfgScalar densityMax = *std::max_element(densities.begin(), densities.end());
    cfgScalar densityMin = *std::min_element(densities.begin(), densities.end());
    std::cout << "density min = " << densityMin << "  density max = " << densityMax << std::endl;

    scores = computePropertiesScores(iPoints);
    int ipoint, npoint=(int)scores.size();
    for (ipoint=0; ipoint<npoint; ipoint++)
    {
      //densities[ipoint] = 1-densities[ipoint];
      densities[ipoint] = 1/densities[ipoint];
    }

    //scores = mult(scores, densities);
    scores = densities;

    /*cfgScalar scoreMax = *std::max_element(scores.begin(), scores.end());
    cfgScalar scoreMin = *std::min_element(scores.begin(), scores.end());

    cfgScalar minProba = 0.01;

    if (scoreMax-scoreMin > 1.e-6)
    {
      for (ipoint=0; ipoint<npoint; ipoint++)
      {
        scores[ipoint] = (scores[ipoint]-scoreMin)/(scoreMax-scoreMin);
        scores[ipoint] *= densities[ipoint];

        scores[ipoint] = (1-minProba)*scores[ipoint] + minProba;
      }
    }*/ 
  }
  return scores;
}
