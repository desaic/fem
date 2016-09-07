#include "FamilyExtractorImpl.h"

#include <stddef.h>
#include <assert.h>
#include <fstream>
#include <iostream>

#include "DistanceField.h"
#include "MicrostructureSet.h"
#include "Resampler.h"
#include "cfgUtilities.h"

FamilyExtractorImpl::FamilyExtractorImpl()
{
  m_microstructures = NULL;

  m_reducedDim = 2;

  //m_stage = InitStage;
  m_stage = MLSComputationStage;
}

FamilyExtractorImpl::~FamilyExtractorImpl()
{
}

void FamilyExtractorImpl::init()
{
}

bool FamilyExtractorImpl::step()
{
  bool resOk = true;

  if (m_stage == InitStage)
  {
    init();
    m_stage = PairwiseDistancesComputationStage;
  }
  else if (m_stage == PairwiseDistancesComputationStage)
  {
    runPairwiseDistancesComputation();
    m_stage = MLSComputationStage;
  }
  else if (m_stage == MLSComputationStage)
  {
    runMLSComputationStage();
    m_stage = EndStage;
  }

  return resOk;
}

bool FamilyExtractorImpl::run()
{
  bool resOk = true;
  while (m_stage != EndStage && resOk)
  {
    resOk = step();
  }
  return resOk;
}

cfgScalar FamilyExtractorImpl::computePairwiseDistance(const std::vector<int> &iMatAssignment1, const std::vector<int> &iMatAssignment2)
{
  cfgScalar dist = 0;
  
  assert(iMatAssignment1.size()==iMatAssignment2.size());
  
  int nvoxel = (int)iMatAssignment1.size();
  for (int ivoxel=0; ivoxel<nvoxel; ivoxel++)
  {
    dist += abs(iMatAssignment1[ivoxel]-iMatAssignment2[ivoxel]);
  }
  return dist;
}

void FamilyExtractorImpl::computePairwiseDistances(const MicrostructureSet &iMicrostructures, const std::vector<int> &iPointIndices, std::vector<cfgScalar> &oDistances)
{
  oDistances.clear();

  int npoint = (int)iPointIndices.size();
  for (int ipoint=0; ipoint<npoint; ipoint++)
  {
    int indMicrostructure1 = iPointIndices[ipoint];
    const std::vector<int> &matAssignment1 = iMicrostructures.getMatAssignment(indMicrostructure1);
    for (int jpoint=0; jpoint<ipoint; jpoint++)
    {
      int indMicrostructure2 = iPointIndices[jpoint];
      const std::vector<int> &matAssignment2 = iMicrostructures.getMatAssignment(indMicrostructure2);

      cfgScalar dist = computePairwiseDistance(matAssignment1, matAssignment2);
      oDistances.push_back(dist);
    }
  }
}

void FamilyExtractorImpl::runPairwiseDistancesComputation()
{
  // extract boundary points
  assert(m_microstructures);

  const std::vector<cfgScalar> & parameters = m_microstructures->getAllMatParameters();
  int paramdim = m_microstructures->getParameterDim();

  DistanceField distanceField(paramdim);
  std::vector<cfgScalar> distances = distanceField.computeDistances(parameters);

  int nTargetParticules = 1000;
  cfgScalar minRadius = 0.1;
  std::vector<int> newparticules;
  Resampler resampler;
  resampler.resampleBoundary(minRadius, paramdim, parameters, distances, nTargetParticules, newparticules);
  m_microstructureIndices = newparticules;

  // evaluate pairwise distances
  std::vector<cfgScalar> pairwiseDistances;
  computePairwiseDistances(*m_microstructures, newparticules, pairwiseDistances);

  std::string fileName = "..//..//Output//distances.txt";
  cfgUtil::writeVector2File<cfgScalar>(pairwiseDistances, fileName); 

  std::string fileNameIndices = "..//..//Output//microstructureIndices.txt";
  cfgUtil::writeVector2File<int>(m_microstructureIndices, fileNameIndices); 
}

void FamilyExtractorImpl::runMLSComputationStage()
{
  std::string fileName = "..//..//Output//reducedCoordinates.txt";
  std::vector<cfgScalar> x;
  cfgUtil::readVectorFromFile(fileName, x);

  std::string fileNameIndices = "..//..//Output//microstructureIndices.txt";
  std::vector<int> inds;
  cfgUtil::readVectorFromFile(fileNameIndices, inds);

  m_reducedCoordinates = x;
  m_microstructureIndices = inds;
}

void FamilyExtractorImpl::extractFamilies()
{
  runPairwiseDistancesComputation();
}
