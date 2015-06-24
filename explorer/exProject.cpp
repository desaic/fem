#include "exProject.h"

#include "ElementRegGrid.hpp"
#include "cfgDefs.h"
#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;

exProject::exProject(int idim)
{
  m_blockSize = 3;
  Q_ASSERT(idim==2 || idim==3);
  m_dim = idim;
}

exProject::~exProject()
{
}

void exProject::setPickedStructure(int iStructureIndex, int iStructureLevel, const std::vector<int> &iBaseMaterials) 
{
  m_pickedStructureIndex = iStructureIndex; 
  m_pickedStructureLevel = iStructureLevel;
  m_pickedStructureBaseMaterials = iBaseMaterials;

  emit pickedStructureModified();
}

bool exProject::readBaseMaterialStructures(int iLevel, std::vector<std::vector<int> > &oBaseMaterialStructures)
{
  int level = iLevel;
 
  std::string fileDirectory = m_FileDirectory;
  std::string materialFile = fileDirectory + "Materials_" + std::to_string(iLevel) + ".txt";

  oBaseMaterialStructures.clear();
  readMaterialCombinations(materialFile, oBaseMaterialStructures);

  return true;
}

bool exProject::getMaterialParameters(int iLevel, std::vector<cfgScalar> &oPhysicalParameters, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<std::vector<int> > &oMaterialAssignmentOneCell)
{
  oPhysicalParameters.clear();

  // read data
  int level = iLevel;
  int blockSize = m_blockSize;

  std::string fileDirectory = m_FileDirectory;
  std::string materialFile = fileDirectory + "Materials_" + std::to_string(level) + ".txt";

  std::string stressStrainFilesDir[2];
  int nlevel = 2;
  if (1)
  //if (iLevel<nlevel)
  {
    stressStrainFilesDir[0] =  fileDirectory + "x_level" + std::to_string(level) + "//";
    stressStrainFilesDir[1] =  fileDirectory + "y_level" + std::to_string(level) + "//";
  }
  else
  {
    stressStrainFilesDir[0] =  fileDirectory + "x" + "//";
    stressStrainFilesDir[1] =  fileDirectory + "y" + "//";
  }
  std::string stressStrainBaseFileName = "StressStrain_" + std::to_string(level) + "_" + std::to_string(blockSize);
  int nbFiles = 0;
  if (level==1)
  {
    if (m_dim==3)
    {
      nbFiles = 256;
    }
    else
    {
      nbFiles = 16;
    }
  }
  else if (level==2)
  {
    //nbFiles = 430;
    nbFiles = 65536;
    //nbFiles = 20; 
  }
  else if (level==3)
  {
    nbFiles = 1198;
  }
  std::vector<std::vector<int> > MaterialAssignments, BaseMaterials;
  bool ResOk = computeMaterialParameters(m_dim, materialFile, stressStrainFilesDir, stressStrainBaseFileName, nbFiles, oPhysicalParameters, oBaseMaterialStructures, MaterialAssignments, oMaterialAssignmentOneCell);

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

  std::string fileDirectory = m_FileDirectory;
  std::string stressStrainFilesDir[2];
  stressStrainFilesDir[0] =  fileDirectory + "x_level" + std::to_string(level) + "//";
  stressStrainFilesDir[1] =  fileDirectory + "y_level" + std::to_string(level) + "//";
  std::string stressStrainBaseFileName = "StressStrain_" + std::to_string(level) + "_" + std::to_string(blockSize);

  std::string sampleFile = stressStrainFilesDir[0] + stressStrainBaseFileName + "_" + std::to_string(indComb) + ".txt"; 
  Vector3f forceAxis;
  std::vector<int> materialAssignment;
  std::vector<cfgScalar> stresses, strains;
  ResOk = readData(sampleFile, forceAxis, materialAssignment, stresses, strains);

  int c = 2*blockSize;
  int n[3] = {c*N[0], c*N[1], c*N[2]};

  std::vector<int> cellMaterials;
  if (m_dim==2)
  {
    getMaterialAssignment(n[0], n[1], materialAssignment, N[0], N[1], baseMaterialStructures, cellMaterials);
    n[2] = 1;
  }
  else
  {
    getMaterialAssignment(n[0], n[1], n[2], materialAssignment, N[0], N[1], N[2], baseMaterialStructures, cellMaterials);
  }
  ElementRegGrid * em = new ElementRegGrid(n[0], n[1], n[2]);
  em->me = cellMaterials;

  QSharedPointer<ElementMesh> elementMesh(em);
  return elementMesh;
}
