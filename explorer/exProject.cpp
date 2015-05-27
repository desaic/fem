#include "exProject.h"

#include "ElementRegGrid.hpp"
#include "cfgDefs.h"
#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;

exProject::exProject()
{
  m_blockSize = 3;
}

exProject::~exProject()
{
}

void exProject::setPickedStructure(int iStructureIndex, int iStructureLevel) 
{
  m_pickedStructureIndex = iStructureIndex; 
  m_pickedStructureLevel = iStructureLevel;

  emit pickedStructureModified();
}

bool exProject::readBaseMaterialStructures(int iLevel, std::vector<std::vector<int> > &oBaseMaterialStructures)
{
  int level = iLevel;
 
  std::string fileDirectory = "C://Users//melina//Documents//Projects//fem//vc//runtime//Output//";
  std::string materialFile = fileDirectory + "Materials_" + std::to_string(iLevel) + ".txt";

  oBaseMaterialStructures.clear();
  readMaterialCombinations(materialFile, oBaseMaterialStructures);

  return true;
}

bool exProject::getMaterialParameters(std::vector<Scalar> &oPhysicalParameters, int iLevel)
{
  oPhysicalParameters.clear();

  // read data
  int level = iLevel;
  int blockSize = m_blockSize;

  std::string fileDirectory = "C://Users//melina//Documents//Projects//fem//vc//runtime//Output//";
  std::string materialFile = fileDirectory + "Materials_" + std::to_string(level) + ".txt";

  std::string stressStrainFilesDir[2];
  stressStrainFilesDir[0] =  fileDirectory + "x_level" + std::to_string(level) + "//";
  stressStrainFilesDir[1] =  fileDirectory + "y_level" + std::to_string(level) + "//";

  std::string stressStrainBaseFileName = "StressStrain_" + std::to_string(level) + "_" + std::to_string(blockSize);
  int nbFiles = 0;
  if (level==1)
  {
    nbFiles = 256;
  }
  else if (level==2)
  {
    nbFiles = 430;
  }
  bool ResOk = computeMaterialParameters(materialFile, stressStrainFilesDir, stressStrainBaseFileName, nbFiles, oPhysicalParameters);

  return true;
}

QSharedPointer<ElementMesh> exProject::computeElementMesh(int iCombIndex, int iLevel)
{
  std::vector<std::vector<int> > baseMaterialStructures;
  bool ResOk = readBaseMaterialStructures(iLevel,baseMaterialStructures);

  int side = pow(2,iLevel-1);
  int N[3] = {side,side,side};

  // level 1
  int level = iLevel;
  int blockSize = m_blockSize;
  int indComb = iCombIndex;

  std::string fileDirectory = "C://Users//melina//Documents//Projects//fem//vc//runtime//Output//";
  std::string stressStrainFilesDir[2];
  stressStrainFilesDir[0] =  fileDirectory + "x_level" + std::to_string(level) + "//";
  stressStrainFilesDir[1] =  fileDirectory + "y_level" + std::to_string(level) + "//";
  std::string stressStrainBaseFileName = "StressStrain_" + std::to_string(level) + "_" + std::to_string(blockSize);

  std::string sampleFile = stressStrainFilesDir[0] + stressStrainBaseFileName + "_" + std::to_string(indComb) + ".txt"; 
  Vector3f forceAxis;
  std::vector<int> materialAssignment;
  std::vector<Scalar> stresses, strains;
  ResOk = readData(sampleFile, forceAxis, materialAssignment, stresses, strains);

  int c = 2*blockSize;
  int n[3] = {c*N[0], c*N[1], c*N[2]};

  std::vector<int> cellMaterials;
  getMaterialAssignment(n[0], n[1], n[2], materialAssignment, N[0], N[1], N[2], baseMaterialStructures, cellMaterials);

  ElementRegGrid * em = new ElementRegGrid(n[0], n[1], n[2]);
  em->me = cellMaterials;

  QSharedPointer<ElementMesh> elementMesh(em);
  return elementMesh;
}
