#ifndef FamilyExtractorImpl_h
#define FamilyExtractorImpl_h

#include "FamilyExtractor.h"
#include "cfgDefs.h"

class MicrostructureSet;

class FamilyExtractorImpl: public FamilyExtractor
{
public:
  FamilyExtractorImpl();
  virtual ~FamilyExtractorImpl();

  virtual void setOption(int iOption){m_option = iOption;};
  virtual void setMicrostructures(MicrostructureSet *iMicrostructures, const std::vector<int> *iIndices=NULL) {m_microstructures = iMicrostructures; m_indices = *iIndices;}

  virtual bool step();
  virtual bool run();

  virtual const std::vector<int> & getMicrostructureIndices() {return m_microstructureIndices;}
  virtual const std::vector<cfgScalar> & getReducedCoordinates() {return m_reducedCoordinates;}
  virtual int getReducedDimension() {return m_reducedDim;}

private:
  enum Stage
  {
    InitStage,
    PairwiseDistancesComputationStage,
    MLSComputationStage,
    EndStage
  };

private:
  void init();
  void extractFamilies();
  void runPairwiseDistancesComputation();
  void runMLSComputationStage();

  cfgScalar computePairwiseDistance(const std::vector<int> &iMatAssignment1, const std::vector<int> &iMatAssignment2);
  void computePairwiseDistances(const MicrostructureSet &iMicrostructures, const std::vector<int> &iPointIndices, std::vector<cfgScalar> &oDistances);

  // Utilities
  bool writeLowerLeftTriangularMatrix(const std::string &iFileName, std::vector<cfgScalar> &iValues, int iMatrixDim, bool iHasValuesOnDiag);

private:
  Stage m_stage;

  MicrostructureSet * m_microstructures;
  std::vector<int> m_indices;

  std::vector<int> m_microstructureIndices;
  std::vector<cfgScalar> m_reducedCoordinates;
  int m_reducedDim;

  int m_option;
};

#endif