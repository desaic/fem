#ifndef SamplesGeneratorImpl_h
#define SamplesGeneratorImpl_h

#include "cfgDefs.h"

class MaterialQuad;
class MaterialQuad2D;
class ElementRegGrid;
class ElementRegGrid2D;
class Stepper;
class Stepper2D;
class Vector3f;

class SamplesGeneratorImpl
{
public:
  SamplesGeneratorImpl(int iDim);
  ~SamplesGeneratorImpl();

  void setOutputDirectory(const std::string &iDirectory) {m_OutputDirectory = iDirectory;}
  const std::string & getOutputDirectory() {return m_OutputDirectory;}

  int run();

private:
  int computeMaterialParameters(std::string & iStepperType, int iLevel, int iNbCombPrevLevel, bool iReadSingleFile, bool iWriteSingleFile, bool iReadFullDeformation, bool iWriteFullDeformation);
  int computeMaterialParametersLevel1(std::string & iStepperType, bool iWriteSingleFile, bool iWriteFullDeformation);
  int computeMaterialParameters(std::string iStepperType,  const std::vector<std::vector<int> > &iMaterialAssignments, 
                                const std::vector<std::vector<int> > &iBaseMaterialStructures, int iBaseMatStructureSize[3], int iLevel, int iBlockRep, bool iWriteSingleFile, bool iWriteFullDeformation);

  void vectorMat2MatrixMat(const std::vector<int> &iVecMaterials, int oMatrixMaterials[2][2][2]);
  void assignMaterials(ElementRegGrid * iElementGrid, int iMacroStructureMaterials[2][2][2]);

  void int2CombRandom(std::vector<int> &oMaterials, int nMat, int nVox);
  void int2Comb_biMaterialStructures(std::vector<int> &oMaterials, int idx, int nMat, int iN[3]);
  void int2Comb(std::vector<int> &oMaterials, int idx, int nMat, int nVox);
  bool sampleDeformation(int iN[3], std::vector<MaterialQuad> &iMaterial, const std::string iStepperType, float iMaxForce, int iForceAxis, int iNumberOfSample, const std::vector<int> & iMaterials,
                         std::vector<float> &oStresses, std::vector<float> &oStrains);
  bool sampleDeformation(int iN[2], std::vector<MaterialQuad2D> &iMaterial, const std::string iStepperType, float iMaxForce, int iForceAxis, int iNumberOfSample, const std::vector<int> & iMaterials,
                         std::vector<float> &oStresses, std::vector<float> &oStrains);
  bool sampleDeformation(int iN[2], std::vector<MaterialQuad2D> &iMaterial, const std::string iStepperType, float iMaxForce, int iForceAxis, int iNumberOfSample, const std::vector<int> & iMaterials,
                         std::vector<float> &oStresses, std::vector<std::vector<float> > &oDeformations);

  void sampleMaterialSpace();
  int checkBorderSimilarity(const std::vector<int> &iMaterialAssignment, int N[3], const std::vector<std::vector<int> > &iBaseMaterialStructures);

  bool computeForceDeformationSample(ElementRegGrid * iElementGrid, std::string iStepperType, int iNumberOfSteps, bool &oConverged);
  bool computeForceDeformationSample(ElementRegGrid2D * iElementGrid, std::string iStepperType, int iNumberOfSteps, bool &oConverged);
  bool computeEquilibrium(Stepper * iStepper, ElementRegGrid * iElementGrid, int iNumberOfSteps, bool &oConverged);
  bool computeEquilibrium(Stepper2D * iStepper, ElementRegGrid2D * iElementGrid, int iNumberOfSteps, bool &oConverged);
  ElementRegGrid * createPhysicalSystem(int iN[3], Vector3f iForce, std::vector<MaterialQuad> &iMaterials);
  ElementRegGrid * createPhysicalSystem(int iN[3], std::vector<MaterialQuad> &iMaterials);
  ElementRegGrid2D * createPhysicalSystem(int iN[2], std::vector<MaterialQuad2D> &iMaterials);
  void setExternalForces(ElementRegGrid * iElementGrid, int iAxis, int iSide, float iForceMagnitude);
  void setExternalForces(ElementRegGrid2D * iElementGrid, int iAxis, int iSide, float iForceMagnitude);
  float computeStrain(const ElementRegGrid * iElementGrid, int iAxis);
  float computeStrain(const ElementRegGrid2D * iElementGrid, int iAxis);

  void getNewCombinations(const std::vector<int> &iBoundaryPoints, const std::vector<cfgScalar> &iPoints, int N[3], const std::vector<std::vector<int> > &iBaseMaterialStructures, int iNbCombinations,
                          std::vector<std::vector<int> > &oNewMaterialAssignments, std::vector<int> &oNewBaseMaterialStructures);
  void getNewCombinationsV2(const std::vector<int> &iBoundaryPoints, const std::vector<cfgScalar> &iPoints, int N[3], int iNbCombinations, std::vector<std::vector<int> > &oNewMaterialAssignments, std::vector<int> &oNewBaseMaterialStructures);
  void getNewCombinationsV3(const std::vector<int> &iBoundaryPoints, const std::vector<cfgScalar> &iPoints, int N[3], int iNbCombinations, std::vector<std::vector<int> > &oNewMaterialAssignments, std::vector<int> &oNewBaseMaterialStructures);
  void getNewCombinationsV4(const std::vector<int> &iBoundaryPoints, const std::vector<cfgScalar> &iPoints, int N[3], int iNbCombinations, std::vector<std::vector<int> > &oNewMaterialAssignments, 
                            const std::vector<std::vector<int> >  &iNewBaseMaterialStructures);
  void getNewCombinationsV5(const std::vector<int> &iBoundaryPoints, const std::vector<cfgScalar> &iPoints, int N[3], int iNbCombinations, std::vector<std::vector<int> > &oNewMaterialAssignments, 
                            const std::vector<std::vector<int> >  &iNewBaseMaterialStructures);

  int computeMaterialParametersIncremental(std::string & iStepperType, int iLevel);
  int computeMaterialParameters(std::string iStepperType , const std::vector<std::vector<int> > &iMaterialAssignments, int iNewMatStructureSize[3],
                                const std::vector<std::vector<int> > &iBaseMaterialStructures, int iBaseMatStructureSize[3], int iLevel, int iBlockRep, bool iWriteSingleFile, bool iWriteFullDeformation, const std::string &iPostfix ="");

  void growStructure(int N[2], const std::vector<std::vector<int> > &materialAssignments, std::vector<std::vector<int> > &oNewMaterialAssignments);

  int getClosestPoint(Vector3f &iP, const std::vector<cfgScalar> &iPoints);
  int getRandomMatWithLinearProba(std::multimap<float, int> &iDist2MatIndex, float eps, float r);

  int computeVariations(std::string & iStepperType, int iLevel);

private:
  std::string m_OutputDirectory;
  bool m_UseLinearMaterial;
  int m_dim;
  int m_blockRep;
  int m_nbSubdivisions;
  std::vector<MaterialQuad> m_mat;
  std::vector<MaterialQuad2D> m_mat2D;
};

#endif





