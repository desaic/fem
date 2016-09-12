#include "FamilyExtractorImpl.h"

#include <stddef.h>
#include <assert.h>
#include <fstream>
#include <iostream>

#include "MicrostructureSet.h"
#include "cfgUtilities.h"
#include "cfgMaterialUtilities.h"

FamilyExtractorImpl::FamilyExtractorImpl()
{
  m_microstructures = NULL;

  m_reducedDim = 2;

  m_stage = InitStage;
  //m_stage = MLSComputationStage;
}

FamilyExtractorImpl::~FamilyExtractorImpl()
{
}

void FamilyExtractorImpl::init()
{
  if (m_microstructures && m_indices.size()==0) 
  {
    int nstructure = m_microstructures->getNumStructures();
    m_indices = cfgMaterialUtilities::genIncrementalSequence(0, nstructure-1, 1);
  }
}

bool FamilyExtractorImpl::step()
{
  bool resOk = true;

  if (m_stage == InitStage)
  {
    init();
    if (m_option==0)
    {
      m_stage = PairwiseDistancesComputationStage;
    }
    else
    {
      m_stage = MLSComputationStage;
    }
  }
  else if (m_stage == PairwiseDistancesComputationStage)
  {
    runPairwiseDistancesComputation();
    m_stage = EndStage;
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

  // evaluate pairwise distances
  std::vector<cfgScalar> pairwiseDistances;
  computePairwiseDistances(*m_microstructures, m_indices, pairwiseDistances);

  std::string fileName = "..//..//Output//distances.txt";
  cfgUtil::writeVector2File<cfgScalar>(pairwiseDistances, fileName); 

  std::string fileNameIndices = "..//..//Output//microstructureIndices.txt";
  cfgUtil::writeVector2File<int>(m_indices, fileNameIndices); 

  m_microstructureIndices = m_indices;
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
