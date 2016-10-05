#ifndef FamilyGeneratorImpl_h
#define FamilyGeneratorImpl_h

#include "FamilyGenerator.h"
#include "cfgDefs.h"

class FamilyGeneratorImpl: public FamilyGenerator
{
public:
  FamilyGeneratorImpl();
  virtual ~FamilyGeneratorImpl();

  virtual void setMicrostructureSize(int iDim, int n[3]);

  virtual bool step();
  virtual bool run();

  virtual const std::vector<std::vector<int> > & getMicrostructures() {return m_matAssignments;}
  virtual void getSize(int oN[3]) {for (int i=0; i<3; i++) {oN[i] = m_n[i];}}

private:
  enum Stage
  {
    InitStage,
    MicrostructureGenerationStage,
  //  PairwiseDistancesComputationStage,
  //  MLSComputationStage,
    EndStage
  };

private:
  void init();

  void genMicrostructure_Type1(const cfgScalar &iSquareLength, const Vector2S &iRec1_Lengths, const Vector2S &iRec2_Lengths, const Vector2S &iRec3_Lengths, std::vector<int> &oMicrostructure);

  // x^2/a^2 + y^2/b^2 + z^2/c^2= 1 params = [a b c]
  void addEllipse(cfgScalar *iParams, cfgScalar *iCenterCoords, cfgScalar iThickness, cfgScalar iRotationAngle, int iMaterial, std::vector<int> &ioMatAssignment);

private:
  Stage m_stage;

  int m_dim;
  int m_n[3];

  std::vector<std::vector<int> > m_matAssignments;
};

#endif