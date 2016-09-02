#ifndef FamilyExtractorImpl_h
#define FamilyExtractorImpl_h

#include "FamilyExtractor.h"

class MicrostructureSet;

class FamilyExtractorImpl: public FamilyExtractor
{
public:
  FamilyExtractorImpl();
  virtual ~FamilyExtractorImpl();

  virtual void setMicrostructures(MicrostructureSet *iMicrostructures) {m_microstructures = iMicrostructures;}

  virtual bool step();
  virtual bool run();

private:
  enum Stage
  {
    InitStage,
    ComputationStage,
    EndStage
  };

private:
  void init();
  void extractFamilies();

private:
  Stage m_stage;

  MicrostructureSet * m_microstructures;
};

#endif