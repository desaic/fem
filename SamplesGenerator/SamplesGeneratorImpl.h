#ifndef SamplesGeneratorImpl_h
#define SamplesGeneratorImpl_h

#include "cfgDefs.h"
#include <set>

class MaterialQuad;
class MaterialQuad2D;
class ElementRegGrid;
class ElementRegGrid2D;
class Stepper;
class Stepper2D;

class SamplesGeneratorImpl
{
public:
  SamplesGeneratorImpl();
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
  bool sampleDeformation(int iN[3], std::vector<MaterialQuad> &iMaterial, const std::string iStepperType, float iMaxForce, int iForceAxis, int iNumberOfSample, const std::vector<int> & iMaterials,
                          std::vector<float> &oStresses, std::vector<std::vector<float> > &oDeformations);
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
  ElementRegGrid * createPhysicalSystem(int iN[3], Vector3S iForce, std::vector<MaterialQuad> &iMaterials);
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

  int computeMaterialParametersIncremental(std::string & iStepperType, int iLevel, std::vector<std::vector<int> > &oNewMaterialAssignments);
  int computeMaterialParameters(std::string iStepperType , const std::vector<std::vector<int> > &iMaterialAssignments, int iNewMatStructureSize[3],
                                const std::vector<std::vector<int> > &iBaseMaterialStructures, int iBaseMatStructureSize[3], int iLevel, int iBlockRep, bool iWriteSingleFile, bool iWriteFullDeformation, const std::string &iPostfix ="");
  int computeDeformation(std::string iStepperType , const std::vector<std::vector<int> > &iMaterialAssignments, int iMatStructureSize[3], const std::vector<std::vector<int> > &iBaseMaterialStructures, int iBaseMatStructureSize[3],
                        std::vector<std::vector<float> > oStresses[3], std::vector<std::vector<std::vector<float> > > oX[3]);

  void growStructure(int N[2], const std::vector<std::vector<int> > &materialAssignments, std::vector<std::vector<int> > &oNewMaterialAssignments);
  void growStructureDoubleSize(int N[2], const std::vector<std::vector<int> > &materialAssignments, std::vector<std::vector<int> > &oNewMaterialAssignments, int iRowToDuplicate=-1, int iColToDuplicate=-1);
  void growStructureDoubleSize(int N[2], const std::vector<int> &matAssignment, std::set<std::vector<int> > &ioNewMaterialAssignments, int iRowToDuplicate, int iColToDuplicate);
  void growStructureDoubleSize(int N[3], const std::vector<int> &matAssignment, std::set<std::vector<int> > &ioNewMaterialAssignments, int iLayerXToDuplicate, int iLayerYToDuplicate, int iLayerZToDuplicate);

  int getClosestPoint(Vector3S &iP, const std::vector<cfgScalar> &iPoints);
  int getRandomMatWithLinearProba(std::multimap<float, int> &iDist2MatIndex, float eps, float r);

  int computeVariations(std::string & iStepperType, int iLevel);
  int computeVariationsV2(int iLevel, int iSubLevel, const std::vector<std::vector<std::vector<int> > > &iMaterialAssignmentsPreviousLevels, const std::vector<std::vector<std::vector<int> > > &iBaseMaterialStructuresPreviousLevels, 
                          const std::vector<std::vector<float> > &iPhysicalParametersPreviousLevels, std::vector<std::vector<int> > &oNewMaterialAssignments);
  int computeVariationsV3(int iLevel, int iSubLevel, const std::vector<std::vector<std::vector<int> > > &iMaterialAssignmentsPreviousLevels, const std::vector<std::vector<std::vector<int> > > &iBaseMaterialStructuresPreviousLevels, 
                          const std::vector<std::vector<float> > &iPhysicalParametersPreviousLevels, std::vector<std::vector<int> > &oNewMaterialAssignments);
  void computeVariationsSMC(int iLevel, std::string & iStepperType, const std::vector<std::vector<int> > &iMaterialAssignments, const std::vector<std::vector<int> > &iBaseMaterialStructures, 
                           const std::vector<float> &iPhysicalParametersPreviousLevels, const std::vector<float> &iTensorsPreviousLevels, bool iGrowStructures,
                           std::vector<std::vector<int> > &oNewMaterialAssignments, std::vector<cfgScalar> &oNewPhysicalParameters, std::vector<cfgScalar> &oNewTensors);

  void writeStressDeformationFile(const std::vector<std::vector<int> > &iMaterialAssignments, int iMatStructureSize[3], const std::vector<std::vector<float> > iStresses[2], const std::vector<std::vector<std::vector<float> > > iX[2], std::string iPostfix="");

  void computeParameters(std::string & iStepperType, const std::vector<std::vector<int> > &iMaterialAssignments, int n[3], const std::vector<std::vector<int> > &iBaseMaterialAssignments, int N[3], 
                         std::vector<cfgScalar> &oParameters, std::vector<std::vector<std::vector<float> > > oX[3]);
  void writeFiles(int iLevel, const std::vector<std::vector<int> > &iMaterialAssignments, const std::vector<std::vector<int> > &iBaseMaterialStructures, const std::vector<cfgScalar> &iParameters, const std::vector<std::vector<std::vector<float> > > iX[2],
                  const std::string iPostFix="");
  void writeFiles(int iLevel, const std::vector<std::vector<int> > &iMaterialAssignments, const std::vector<std::vector<int> > &iBaseMaterialStructures, const std::vector<cfgScalar> &iParameters,
                  const std::vector<cfgScalar> &iTensorsValues, const std::string iPostFix="");
  void writeFiles(int iLevel, const std::vector<std::vector<float> > &iMaterialDistributions, const std::vector<cfgScalar> &iParameters, const std::vector<cfgScalar> &iTensorsValues, const std::string iPostFix="");

  bool readFiles(int iLevel, std::vector<std::vector<int> > &oMaterialAssignments, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<cfgScalar> &oParameters);
  bool readFiles(int iLevel, std::vector<std::vector<int> > &oMaterialAssignments, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<cfgScalar> &oParameters, std::vector<cfgScalar> &oTensors, const std::string iPostFix="");
  void concatenateFiles(int iLevel, int indexMin, int indexMax, const std::string iPostFix="", const std::string iNewPostFix="");
  void concatenateFiles(int iLevel, const std::string iPostFix1, const std::string iPostFix2, const std::string iNewPostFix="");

  int getNbParameters();
  int getNbFreeCells(int n[3]);
  std::vector<int> getFreeCellsMaterialAssignment(int n[3], const std::vector<int> &iMatAssignments);
  std::vector<int>  getFullMaterialAssignment(int n[3], const std::vector<int> &iFreeCellMatAssignments);
  void computeParametersAndTensorValues(int n[3], const std::vector<std::vector<int> > &iMatAssignments, std::vector<cfgScalar> &oParameters, std::vector<cfgScalar> &oTensorValues);
  void computeParametersAndTensorValues(int n[3], const std::vector<std::vector<double> > &iMatDistributions, std::vector<cfgScalar> &oParameters, std::vector<cfgScalar> &oTensorValues);

  void runContinuousOptimization(int iLevel, int iCycle, bool iFixNonManifoldStructure, const std::vector<std::vector<int> > &iMaterialAssignments, 
                                 const std::vector<float> &iParameters, const std::vector<float> &iTensors, 
                                 std::vector<std::vector<int> > &oNewMaterialAssignments, std::vector<cfgScalar> &oNewPhysicalParameters, std::vector<cfgScalar> &oNewTensors);

  void runContinuousOptimization(int iLevel, int iCycle, bool iFixNonManifoldStructure, const std::vector<std::vector<int> > &iMaterialAssignments, const std::vector<float> &iParameters,
                                 const std::vector<float> &iTensors, const std::vector<float> &iOffsetedParameters, const std::vector<int> &iNewParticules, std::vector<std::vector<int> > &oNewMaterialAssignments,
                                 std::vector<cfgScalar> &oNewPhysicalParameters, std::vector<cfgScalar> &oNewTensors);

  bool fixNonManifoldStructure(int n[3], std::vector<int> &ioMatAssignments);
  bool fixNonManifoldStructure2D(int n[2], std::vector<int> &ioMatAssignments);

  bool readDisneyFiles(const std::string &iFileName, int iDim, bool iCubic, bool iOrthotropic, std::vector<cfgScalar> &ioPhysicalParameters);
  bool readDisneyFiles(int iDim, bool iCubic, bool iOrthotropic, std::vector<cfgScalar> &oPhysicalParameters);

private:
  std::string m_OutputDirectory;
  bool m_UseLinearMaterial;
  int m_dim;
  int m_blockRep;
  int m_nbSubdivisions;
  std::vector<MaterialQuad> m_mat;
  std::vector<MaterialQuad2D> m_mat2D;
  bool m_orthotropicOnly;
  bool m_cubicOnly;
  bool m_continuousMatDist;
  bool m_filterOutDisconnetedStructures;
};

#endif





