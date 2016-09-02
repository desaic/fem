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

void Resampler::resample(cfgScalar iMinRadius, int iDim, const std::vector<cfgScalar> &iPoints, const std::vector<cfgScalar> &iScores, int iNbTargetParticules, std::vector<int> &oParticules, cfgScalar iRatio)
{
  oParticules.clear();

  int ntargetparticules = iNbTargetParticules;
  cfgScalar ratio = iRatio;
  if (ratio>0)
  {
    int ntop = ratio*iNbTargetParticules;
    std::multimap<cfgScalar, int> score2Particule;
    int npoint = (int)iScores.size();
    for (int ipoint=0; ipoint<npoint; ipoint++)
    {
      score2Particule.insert(std::make_pair(iScores[ipoint], ipoint));
    }
    int ind=0;
    std::multimap<cfgScalar,int>::reverse_iterator it;
    for (it=score2Particule.rbegin(); it!=score2Particule.rend() && ind<npoint; it++, ind++)
    {
      oParticules.push_back(it->second);
    }
    ntargetparticules -= ntop;
  }
  std::vector<int> pointIndices;
  if (iMinRadius>0)
  {
    resampleData(iMinRadius, iDim, iPoints, iScores, pointIndices);
  }
  else
  {
    pointIndices = genIncrementalSequence(0, (int)iScores.size()-1);
  }
   // Systematic resampling;
  std::vector<double> scores = convertVec<cfgScalar, double>(getSubVector(iScores, pointIndices));

  std::vector<double> partialSums(scores.size());
  std::partial_sum(scores.begin(), scores.end(), &partialSums[0]);

  partialSums = div<double>(partialSums, partialSums.back());

  double r0 = (double)rand()/(cfgScalar)RAND_MAX;
  int indParticule = 0;

  for (int iparticule=0; iparticule<ntargetparticules; iparticule++)
  {
    double r = (ntargetparticules*r0 + iparticule)/(double)(ntargetparticules);
    if (r>1)
    {
      r -= 1;
      indParticule = 0;
    }
    double upperBound = partialSums[indParticule];
    while (r > upperBound)
    {
      indParticule++;
      upperBound = partialSums[indParticule];
    }
    oParticules.push_back(pointIndices[indParticule]);
  }
}

void Resampler::resample(const std::vector<cfgScalar> &iScores, int iNbTargetParticules, std::vector<int> &oParticules, cfgScalar iRatio)
{
  oParticules.clear();

  int ntargetparticules = iNbTargetParticules;
  cfgScalar ratio = iRatio;
  if (ratio>0)
  {
    int ntop = ratio*iNbTargetParticules;
    std::multimap<cfgScalar, int> score2Particule;
    int npoint = (int)iScores.size();
    for (int ipoint=0; ipoint<npoint; ipoint++)
    {
      score2Particule.insert(std::make_pair(iScores[ipoint], ipoint));
    }
    int ind=0;
    std::multimap<cfgScalar,int>::reverse_iterator it;
    for (it=score2Particule.rbegin(); it!=score2Particule.rend() && ind<npoint; it++, ind++)
    {
      oParticules.push_back(it->second);
    }
    ntargetparticules -= ntop;
  }


  std::vector<double> scores = convertVec<float, double>(iScores);

  // Systematic resampling;
  std::vector<double> partialSums(scores.size());
  std::partial_sum(scores.begin(), scores.end(), &partialSums[0]);

  partialSums = div<double>(partialSums, partialSums.back());

  double r0 = (double)rand()/(cfgScalar)RAND_MAX;
  int indParticule = 0;

  for (int iparticule=0; iparticule<ntargetparticules; iparticule++)
  {
    double r = (ntargetparticules*r0 + iparticule)/(double)(ntargetparticules);
    if (r>1)
    {
      r -= 1;
      indParticule = 0;
    }
    double upperBound = partialSums[indParticule];
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
  while (particules.size()<ntargetparticules)
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

void Resampler::resampleUsingNormalDistribution(const std::vector<cfgScalar> &iDistances, int iNbTargetParticules, std::vector<int> &oParticules)
{
  oParticules.clear();

  cfgScalar min = getNegativeMin<cfgScalar>(iDistances);

  cfgScalar l = fabs(min);

  std::vector<cfgScalar> x(iDistances.size());
  int nval = (int)iDistances.size();
  for (int ival=0; ival<nval; ival++)
  {
    cfgScalar dist = iDistances[ival];
    if (dist>0)
    {
      x[ival] = 0;
    }
    else
    {
      x[ival] = fabs(dist);
    }
  }
  std::cout << "max dist = " << l << std::endl;

  std::vector<cfgScalar> scores(nval);
  cfgScalar sigma = 0.05;
  cfgScalar mu = 0;
  for (int ival=0; ival<nval; ival++)
  {
    cfgScalar val = x[ival];
    cfgScalar r = (val-mu)/sigma;
    cfgScalar p = 1./(sigma*sqrt(2*M_PI)) * exp(-0.5*r*r);
    scores[ival] = p;
  }
  resample(scores, iNbTargetParticules, oParticules);
}


