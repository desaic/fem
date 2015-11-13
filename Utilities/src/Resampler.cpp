/*=====================================================================================*/
/*! \file		Resampler.cpp
	\author		skourasm
	\brief		Implementation of class Resampler
 */
/*=====================================================================================*/

#include "Resampler.h"

#include <cfgUtilities.h>
using namespace cfgUtil;

Resampler::Resampler()
{
}

Resampler::~Resampler()
{
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


