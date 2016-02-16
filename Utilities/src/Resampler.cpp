/*=====================================================================================*/
/*! \file		Resampler.cpp
	\author		skourasm
	\brief		Implementation of class Resampler
 */
/*=====================================================================================*/

#include "Resampler.h"

#include "DistanceTool.h"

#include <cfgUtilities.h>
using namespace cfgUtil;

#include <cfgMaterialUtilities.h>
using namespace cfgMaterialUtilities;

Resampler::Resampler()
{
}

Resampler::~Resampler()
{
}

void Resampler::resample(cfgScalar iMinRadius, int iDim, const std::vector<cfgScalar> &iPoints, const std::vector<cfgScalar> &iScores, int iNbTargetParticules, std::vector<int> &oParticules)
{
  oParticules.clear();

  std::vector<int> pointIndices;
  if (iMinRadius>0)
  {
    resampleData(iMinRadius, iDim, iPoints, iScores, pointIndices);
  }
  else
  {
    pointIndices = genIncrementalSequence(0, (int)iScores.size()-1);
  }
  std::vector<cfgScalar> scores = getSubVector(iScores, pointIndices);
  normalize(scores);

  // Systematic resampling;
  std::vector<cfgScalar> partialSums(scores.size());
  std::partial_sum(scores.begin(), scores.end(), &partialSums[0]);

  cfgScalar r0 = (cfgScalar)rand()/(cfgScalar)RAND_MAX;
  cfgScalar step = (cfgScalar)1/(cfgScalar)(iNbTargetParticules);
  int indParticule = 0;

  for (int iparticule=0; iparticule<iNbTargetParticules; iparticule++)
  {
    cfgScalar r = r0 + iparticule*step;
    if (r>1)
    {
      r -= 1;
      indParticule = 0;
    }
    cfgScalar upperBound = partialSums[indParticule];
    while (r > upperBound)
    {
      indParticule++;
      upperBound = partialSums[indParticule];
    }
    oParticules.push_back(pointIndices[indParticule]);
  }
}

void Resampler::resample(const std::vector<cfgScalar> &iScores, int iNbTargetParticules, std::vector<int> &oParticules)
{
  oParticules.clear();


  std::vector<cfgScalar> scores = iScores;
  normalize(scores);

  // Systematic resampling;
  std::vector<cfgScalar> partialSums(scores.size());
  std::partial_sum(scores.begin(), scores.end(), &partialSums[0]);

  cfgScalar r0 = (cfgScalar)rand()/(cfgScalar)RAND_MAX;
  cfgScalar step = (cfgScalar)1/(cfgScalar)(iNbTargetParticules);
  int indParticule = 0;

  for (int iparticule=0; iparticule<iNbTargetParticules; iparticule++)
  {
    cfgScalar r = r0 + iparticule*step;
    if (r>1)
    {
      r -= 1;
      indParticule = 0;
    }
    cfgScalar upperBound = partialSums[indParticule];
    while (r > upperBound)
    {
      indParticule++;
      upperBound = partialSums[indParticule];
    }
    oParticules.push_back(indParticule);
  }


  /*int nInitParticules = scores.size();
  int indParticule = 0;

  std::vector<int> particules;
  while (particules.size()<iNbTargetParticules)
  {
    float r = (float)rand()/float(RAND_MAX);
    int indParticule = rand() % nInitParticules;
    float p = iScores[indParticule];
    if (p >= r)
    {
      particules.push_back(indParticule);
    }
    //indParticule  = indParticule % nInitParticules;
    //particules.push_back(indParticule++);
  }
  oParticules = particules; */ 
}


void Resampler::normalize(std::vector<cfgScalar> &ioValues)
{
  cfgScalar sumValues = std::accumulate(ioValues.begin(), ioValues.end(), (cfgScalar)0);
  ioValues = mult(ioValues, 1/sumValues);
}

void Resampler::resampleBoundary(cfgScalar iMinRadius, int iDim, const std::vector<cfgScalar> &iPoints, const std::vector<cfgScalar> &iScores, int iNbTargetParticules, std::vector<int> &oParticules)
{
  oParticules.clear();

  std::set<int> pointsToDiscard;

  std::vector<cfgScalar> initPoints = iPoints;
  std::vector<cfgScalar> lengths(iDim, 1);
  bool resOk = rescaleData(initPoints, iDim, lengths);
  std::vector<double> points = convertVec<cfgScalar, double>(initPoints);

  DistanceTool distanceTool(points, iDim);

  std::vector<int> orderedScoreIndices;
  sortValues(iScores, orderedScoreIndices);

  int npoint = (int)orderedScoreIndices.size();
  for (int ipoint=0; ipoint<npoint && oParticules.size()<iNbTargetParticules; ipoint++)
  {
    int indPoint = orderedScoreIndices[npoint-1-ipoint];
    if (pointsToDiscard.count(indPoint)==0)
    {
      oParticules.push_back(indPoint);
      if (iMinRadius>0)
      {
        std::vector<int> closestPoints;
        distanceTool.getClosestPointIndices(&points[iDim*indPoint], iMinRadius, pointsToDiscard, closestPoints);
        for (int iclosestPoint=0; iclosestPoint<(int)closestPoints.size(); iclosestPoint++)
        {
          pointsToDiscard.insert(closestPoints[iclosestPoint]);
        }
      }
    }
  }
}

void Resampler::resampleData(cfgScalar iMinRadius, int iDim, const std::vector<cfgScalar> &iPoints, const std::vector<cfgScalar> &iScores, std::vector<int> &oParticules)
{
  oParticules.clear();

  std::set<int> pointsToDiscard;

  std::vector<cfgScalar> initPoints = iPoints;
  std::vector<cfgScalar> lengths(iDim, 1);
  bool resOk = rescaleData(initPoints, iDim, lengths);
  std::vector<double> points = convertVec<cfgScalar, double>(initPoints);

  DistanceTool distanceTool(points, iDim);

  std::vector<int> orderedScoreIndices;
  sortValues(iScores, orderedScoreIndices);

  int npoint = (int)orderedScoreIndices.size();
  for (int ipoint=0; ipoint<npoint; ipoint++)
  {
    int indPoint = orderedScoreIndices[npoint-1-ipoint];
    if (pointsToDiscard.count(indPoint)==0)
    {
      std::vector<int> closestPoints;
      distanceTool.getClosestPointIndices(&points[iDim*indPoint], iMinRadius, pointsToDiscard, closestPoints);
      for (int iclosestPoint=0; iclosestPoint<(int)closestPoints.size(); iclosestPoint++)
      {
        pointsToDiscard.insert(closestPoints[iclosestPoint]);
      }
      int randomPoint = rand() % closestPoints.size();
      oParticules.push_back(closestPoints[randomPoint]);
    }
  }
}


