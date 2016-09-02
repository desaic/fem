#include "FamilyExtractorImpl.h"

#include <stddef.h>
#include <assert.h>

#include "DistanceField.h"
#include "MicrostructureSet.h"
#include "Resampler.h"

FamilyExtractorImpl::FamilyExtractorImpl()
{
  m_microstructures = NULL;

  m_stage = InitStage;
}

FamilyExtractorImpl::~FamilyExtractorImpl()
{
}

void FamilyExtractorImpl::init()
{
}

void FamilyExtractorImpl::extractFamilies()
{
  // extract boundary points
  assert(m_microstructures);

  const std::vector<cfgScalar> & parameters = m_microstructures->getAllMatParameters();
  int paramdim = m_microstructures->getParameterDim();

  DistanceField distanceField(paramdim);
  std::vector<cfgScalar> distances = distanceField.computeDistances(parameters);

  int nTargetParticules = 100;
  cfgScalar minRadius = 0.1;
  std::vector<int> newparticules;
  Resampler resampler;
  resampler.resampleBoundary(minRadius, paramdim, parameters, distances, nTargetParticules, newparticules);
}

bool FamilyExtractorImpl::step()
{
  bool resOk = true;

  if (m_stage == InitStage)
  {
    init();
    m_stage = ComputationStage;
  }
  else if (m_stage == ComputationStage)
  {
    extractFamilies();
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


