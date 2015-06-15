#ifndef SamplesGeneratorImpl_h
#define SamplesGeneratorImpl_h

#include "cfgDefs.h"

class MaterialQuad;
class ElementRegGrid;
class Stepper;
class Vector3f;

class SamplesGeneratorImpl
{
public:
  SamplesGeneratorImpl();
  ~SamplesGeneratorImpl();

  void setOutputDirectory(const std::string &iDirectory) {m_OutputDirectory = iDirectory;}
  const std::string & getOutputDirectory() {return m_OutputDirectory;}

  int run();

private:
  int computeMaterialParameters(std::vector<MaterialQuad> & iMaterials, std::string & iStepperType, int iLevel, int iNbCombPrevLevel);
  int computeMaterialParametersLevel1(std::vector<MaterialQuad> & iMaterials, std::string & iStepperType);
  int computeMaterialParameters(std::vector<MaterialQuad> &iMaterials, std::string iStepperType,  const std::vector<std::vector<int> > &iMaterialAssignments, 
                                const std::vector<std::vector<int> > &iBaseMaterialStructures, int iBaseMatStructureSize[3], int iLevel, int iBlockRep);

  void getMaterialAssignment(int nx, int ny, int nz, const std::vector<int> &iMaterialCombIndices, int nX, int nY, int nZ, const std::vector<std::vector<int> > &iBaseCellMaterials, std::vector<int> &oMaterials);
  void vectorMat2MatrixMat(const std::vector<int> &iVecMaterials, int oMatrixMaterials[2][2][2]);
  void assignMaterials(ElementRegGrid * iElementGrid, int iMacroStructureMaterials[2][2][2]);

  void int2CombRandom(std::vector<int> &oMaterials, int nMat, int nVox);
  void int2Comb_biMaterialStructures(std::vector<int> &oMaterials, int idx, int nMat, int iN[3]);
  void int2Comb(std::vector<int> &oMaterials, int idx, int nMat, int nVox);
  bool sampleDeformation(int iN[3], std::vector<MaterialQuad> &iMaterial, const std::string iStepperType, float iMaxForce, int iForceAxis, int iNumberOfSample, const std::vector<int> & iMaterials,
                         std::vector<float> &oStresses, std::vector<float> &oStrains);

  void sampleMaterialSpace();
  int checkBorderSimilarity(const std::vector<int> &iMaterialAssignment, int N[3], const std::vector<std::vector<int> > &iBaseMaterialStructures);

  bool computeForceDeformationSample(ElementRegGrid * iElementGrid, std::string iStepperType, int iNumberOfSteps, bool &oConverged);
  bool computeEquilibrium(Stepper * iStepper, ElementRegGrid * iElementGrid, int iNumberOfSteps, bool &oConverged);
  ElementRegGrid * createPhysicalSystem(int iN[3], Vector3f iForce, std::vector<MaterialQuad> &iMaterials);
  ElementRegGrid * createPhysicalSystem(int iN[3], std::vector<MaterialQuad> &iMaterials);
  void setExternalForces(ElementRegGrid * iElementGrid, int iAxis, int iSide, float iForceMagnitude);
  float computeStrain(const ElementRegGrid * iElementGrid, int iAxis);

  void getNewCombinations(const std::vector<int> &iBoundaryPoints, const std::vector<cfgScalar> &iPoints, int N[3], const std::vector<std::vector<int> > &iBaseMaterialStructures, int iNbCombinations,
                          std::vector<std::vector<int> > &oNewMaterialAssignments, std::vector<int> &oNewBaseMaterialStructures);
  int getClosestPoint(Vector3f &iP, const std::vector<cfgScalar> &iPoints);

private:
  std::string m_OutputDirectory;
  bool m_UseLinearMaterial;
};

#endif





