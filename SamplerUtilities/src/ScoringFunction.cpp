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

cfgScalar ScoringFunction::estimateDensityRadius(const std::vector<cfgScalar> &iPoints, int iDim)
{
  DensityEstimator densityEstimator(iDim);
  densityEstimator.setPoints(iPoints);
  cfgScalar radius = densityEstimator.estimateRadius();
  return radius;
}

std::vector<cfgScalar> ScoringFunction::computeDensities(const std::vector<cfgScalar> &iPoints, int iDim)
{
  DensityEstimator densityEstimator(iDim);
  densityEstimator.setPoints(iPoints);
  densityEstimator.setKernelRadius(m_densityRadius);
  densityEstimator.init();
  std::vector<cfgScalar> densities;
  if (m_cacheDensities)
  {
    densities = m_cachedUnscaledDensities;
    std::vector<cfgScalar> newDensities = densityEstimator.computeDensities(m_newPointIndexStart, (int)iPoints.size()/iDim, densities);
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

std::vector<cfgScalar> ScoringFunction::computePropertiesScores(const std::vector<cfgScalar> &iPoints, int iDim)
{
  std::vector<cfgScalar> scores;

  int ipoint, npoint=(int)iPoints.size()/iDim;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
  /*  cfgScalar score = -(iPoints[iDim*ipoint]+iPoints[iDim*ipoint+1]+iPoints[iDim*ipoint+2]);
    if (iDim==4)
      score -= iPoints[iDim*ipoint+3];
    scores.push_back(score);
    */ 
    // Favors materials with lower density but higher Young modulus
    cfgScalar score = -iPoints[iDim*ipoint]+iPoints[iDim*ipoint+1]+iPoints[iDim*ipoint+2];
    if (iDim==4)
      score += iPoints[iDim*ipoint+3];
    scores.push_back(score);
  }
  return scores;
}

std::vector<cfgScalar> ScoringFunction::computeScores(const std::vector<cfgScalar> &iPoints)
{
  std::vector<cfgScalar> points = iPoints;

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

  int dim = m_dim;
  if (m_paramsToUse.size()>0)
  {
    int nparam = (int)m_paramsToUse.size();

    points.clear();
    int npoint = (int)iPoints.size()/m_dim;
    for (int ipoint=0; ipoint<npoint; ipoint++)
    {
      for (int iparam=0; iparam<nparam; iparam++)
      {
        int indParam = m_paramsToUse[iparam];
        cfgScalar val = iPoints[m_dim*ipoint+indParam];
        points.push_back(val);
      }
    }
    /*dimToScale.clear();
    if (m_dim==4)
    {
      dimToScale.push_back(0);
    }
    else if (m_dim==5)
    {
      dimToScale.push_back(0);
      dimToScale.push_back(1);
    }*/ 
    dim = nparam;
  }

  std::vector<cfgScalar> scores;
  if (m_useDistanceField)
  {
    std::vector<cfgScalar> lengths(dim, 1);
    if (m_useLogScale)
    {
      cfgScalar eps = 1.e-6;
      convertToLogValues(points, dim, dimToScale, eps);
      rescaleData(points, dim, lengths);
    }
    if (m_densityRadius==0)
    {
      m_densityRadius = estimateDensityRadius(points, dim);
    }
    
    std::cout << "computing distances..." << std::endl;
    DistanceField distanceField(dim);
    scores = distanceField.computeDistances(points);

    bool useDensity = true;
    std::vector<cfgScalar> densities;
    if (useDensity)
    {
      std::cout << "computing densities..." << std::endl;
      densities = computeDensities(points, dim);
   
      bool normalizedDensities = true;
      if (normalizedDensities)
      {
        cfgScalar sigma = 0.05;
        cfgScalar mu = 0;
        for (int ival=0; ival<densities.size(); ival++)
        {
          cfgScalar val = densities[ival];
          cfgScalar r = (val-mu)/sigma;
          cfgScalar p = 1./(sigma*sqrt(2*M_PI)) * exp(-0.5*r*r);
          densities[ival] /= p;
        }
      }

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
    std::vector<cfgScalar> lengths(dim, 1);
    rescaleData(points, dim, lengths);

    if (m_densityRadius==0)
    {
      m_densityRadius = estimateDensityRadius(points, dim);
    }
    std::cout << "density radius = " << m_densityRadius << std::endl;
    std::vector<cfgScalar> densities = computeDensities(points, dim);
    cfgScalar densityMax = *std::max_element(densities.begin(), densities.end());
    cfgScalar densityMin = *std::min_element(densities.begin(), densities.end());
    std::cout << "density min = " << densityMin << "  density max = " << densityMax << std::endl;

    scores = computePropertiesScores(points, dim);
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
