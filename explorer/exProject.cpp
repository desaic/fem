#include "exProject.h"

#include "ElementRegGrid.hpp"
#include "cfgDefs.h"
#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;
#include <cfgUtilities.h>
using namespace cfgUtil;

#include <assert.h>

#include <QStringList>

exProject::exProject(int idim)
{
  m_blockSize = 1;
  Q_ASSERT(idim==2 || idim==3);
  m_dim = idim;
  m_readSingleFile = true;
  m_readFullDeformation = true;
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

  std::string fileDirectory = m_FileDirectory;
  std::string materialFile = fileDirectory + "Materials_" + std::to_string(iLevel) + ".txt";

  oBaseMaterialStructures.clear();
  readMaterialCombinations(materialFile, oBaseMaterialStructures);

  return true;
}

bool exProject::loadFile(const QString &iFileName, int &oLevel, QString &oLabel)
{
  bool resOk = true;
  QString str = iFileName.section("level", 1, 1);
  if (!str.isEmpty())
  {
    QString strLevel = str.section('_', 0, 0);
    int level = strLevel.toInt();

    QStringList substrings = iFileName.split('_');
    QString stringEnd = substrings.last();
    QString stringRoot = iFileName.section(stringEnd, 0, 0);

    std::string fileRootName = stringRoot.toStdString();
    std::string fileExtension = ".bin";

    std::vector<cfgScalar> physicalParameters;
    std::vector<std::vector<int> > baseMaterials, materialAssignments;
    resOk = cfgUtil::readBinary<float>(fileRootName + "params" + fileExtension, physicalParameters);
    if (resOk)
      resOk = cfgUtil::readBinary<int>(fileRootName + "baseMat" + fileExtension, baseMaterials);
    if (resOk)
      resOk = cfgUtil::readBinary<int>(fileRootName + "matAssignments" + fileExtension, materialAssignments);

    if (resOk)
    {
      oLevel = level;

      QString str = iFileName.section("level", 0, 0);
      oLabel = stringRoot;
      oLabel.chop(1);
      oLabel.remove(str);

      addLevelData(level, baseMaterials, materialAssignments, physicalParameters);
    }
  }
  return resOk;
}

void exProject::addLevelData(int ilevel, const std::vector<std::vector<int> > &iBaseMaterials, const std::vector<std::vector<int> > &iMaterialAssignments, std::vector<cfgScalar> &iPhysicalParameters)
{
  m_levels.push_back(ilevel);
  m_baseMaterials.push_back(iBaseMaterials);
  m_materialAssignments.push_back(iMaterialAssignments);
  m_physicalParametersPerLevel.push_back(iPhysicalParameters);

  int defaultVisibility = 1;
  m_levelVisibility.push_back(defaultVisibility);

  emit levelsModified();

  //m_baseMaterials2Indices.clear();
  //m_baseMaterials2Indices.resize(nlevel);
}

bool exProject::getMaterialParameters(int iLevel, std::vector<cfgScalar> &oPhysicalParameters, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<std::vector<int> > &oMaterialAssignmentOneCell)
{
  oPhysicalParameters.clear();

  // read data
  int level = iLevel;
  int blockSize = m_blockSize;
  bool readSingleFile = m_readSingleFile;
  bool readFullDeformation = m_readFullDeformation;

  int nlevel = m_levels.size()+1;

  std::cout << "reading material parameters for level " << iLevel << std::endl;

  std::string fileDirectory = m_FileDirectory;
  int indLevel = 2;
  std::string materialFile = fileDirectory + "Materials_" + std::to_string(indLevel) + ".txt";
  bool readBinary = true;
  std::string extension = (readBinary? ".bin": ".ext");


  bool ResOk = true;
  if (readSingleFile)
  {
    if (readFullDeformation)
    {
      std::string baseFileName = "StressDeformation_" + std::to_string(level) + "_" + std::to_string(blockSize);
      std::string stressDeformationFileNames[2] = {fileDirectory + baseFileName + "_x" + extension, fileDirectory + baseFileName + "_y" + extension};
      ResOk = computeMaterialParametersFromDeformations(materialFile, stressDeformationFileNames, readBinary, oPhysicalParameters, oBaseMaterialStructures, oMaterialAssignmentOneCell);
    }
    else
    {
      std::string baseFileName = "StressStrain_" + std::to_string(level) + "_" + std::to_string(blockSize);
      std::string stressStrainFileNames[2] = {fileDirectory + baseFileName + "_x.txt", fileDirectory + baseFileName + "_y.txt"};

      ResOk = computeMaterialParameters(materialFile, stressStrainFileNames, oPhysicalParameters, oBaseMaterialStructures, oMaterialAssignmentOneCell);
    }
  }
  else
  {
    std::string stressStrainFilesDir[2];
    
    //if (1)
    if (iLevel<nlevel-1)
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
      //nbFiles = 10;
      nbFiles = 65536;
      // nbFiles = 15000; 
    }
    else if (level==3)
    {
      nbFiles = 10000; //;075;
      //nbFiles = 1198;
    }
    ResOk = computeMaterialParameters(materialFile, stressStrainFilesDir, stressStrainBaseFileName, nbFiles, oPhysicalParameters, oBaseMaterialStructures, oMaterialAssignmentOneCell);
  }
  return ResOk;
}

void exProject::setMaterialAssignment(int iLevel, const std::vector<std::vector<int> > &iMaterialAssignment)
{
  assert(iLevel>=0 && iLevel<MAX_LEVELS);

  m_materialAssignments[iLevel] = iMaterialAssignment;
}

QSharedPointer<ElementMesh> exProject::computeElementMesh(int iCombIndex, int iLevel)
{
  std::vector<std::vector<int> > baseMaterialStructures;
  bool ResOk = readBaseMaterialStructures(iLevel,baseMaterialStructures);

  int side = pow(2,iLevel-1);
  int N[3] = {side,side,side};

  // level 1
  int level = iLevel;
  int indComb = iCombIndex;

  std::string fileDirectory = m_FileDirectory;

  std::vector<int> materialAssignment;
  if (m_readSingleFile)
  {
    if (m_materialAssignments[iLevel].size()==0)
    {
      std::string fileRootName = getFileDirectory() + "level" + std::to_string(iLevel) + "_";
      std::string fileExtension = ".bin";
      cfgUtil::readBinary<int>(fileRootName + "matAssignments" + fileExtension, m_materialAssignments[iLevel]);

      /*
      Vector3f forceAxis;
      if (m_readFullDeformation)
      {
        std::string baseFileName = "StressDeformation_" + std::to_string(level) + "_" + std::to_string(m_blockSize);
        std::string stressStrainFileName = fileDirectory + baseFileName + "_x.txt";
        std::vector<std::vector<float> > stresses;
        std::vector<std::vector<std::vector<float> > > x;
        int n[3];
        ResOk = readData(stressStrainFileName,  forceAxis, m_materialAssignments[iLevel], stresses, x, n);
      }
      else
      {
        std::string baseFileName = "StressStrain_" + std::to_string(level) + "_" + std::to_string(m_blockSize);
        std::string stressStrainFileName = fileDirectory + baseFileName + "_x.txt";
        std::vector<std::vector<float> > stresses, strains;
        ResOk = readData(stressStrainFileName,  forceAxis, m_materialAssignments[iLevel], stresses, strains);
      }*/ 
    }
    materialAssignment = m_materialAssignments[iLevel][indComb];
  }
  else
  {
    
    std::string stressStrainFilesDir[2];
    stressStrainFilesDir[0] =  fileDirectory + "x_level" + std::to_string(level) + "//";
    stressStrainFilesDir[1] =  fileDirectory + "y_level" + std::to_string(level) + "//";
    std::string stressStrainBaseFileName = "StressStrain_" + std::to_string(level) + "_" + std::to_string(m_blockSize);

    std::string sampleFile = stressStrainFilesDir[0] + stressStrainBaseFileName + "_" + std::to_string(indComb) + ".txt"; 
    Vector3f forceAxis;
    
    std::vector<cfgScalar> stresses, strains;
    ResOk = readData(sampleFile, forceAxis, materialAssignment, stresses, strains);
  }

  int blockSize = 3;
  int c = 2*blockSize;
  int n[3] = {c*N[0], c*N[1], c*N[2]};

  std::vector<int> cellMaterials;
  if (m_dim==2)
  {
    getMaterialAssignment(2, 2, materialAssignment, N[0], N[1], baseMaterialStructures, blockSize, blockSize, 1, cellMaterials);
    n[2] = 1;
  }
  else
  {
    getMaterialAssignment(2, 2, 2, materialAssignment, N[0], N[1], N[2], baseMaterialStructures, blockSize, blockSize, blockSize, 1, cellMaterials);
  }
  ElementRegGrid * em = new ElementRegGrid(n[0], n[1], n[2]);
  em->me = cellMaterials;

  QSharedPointer<ElementMesh> elementMesh(em);
  return elementMesh;
}

QSharedPointer<ElementMesh> exProject::computeElementMeshIncr(int iCombIndex, int iLevel)
{
#if 0
  int N[3] = {iLevel,iLevel,iLevel};

  std::vector<std::vector<int> > baseMaterialStructures;
  bool ResOk = readBaseMaterialStructures(iLevel,baseMaterialStructures);

  int N[3] = {iLevel,iLevel,iLevel};

  // level 1
  int level = iLevel;
  int indComb = iCombIndex;

  std::string fileDirectory = m_FileDirectory;

  std::vector<int> materialAssignment;
  if (m_readSingleFile)
  {
    if (m_materialAssignments[iLevel].size()==0)
    {
      std::string fileRootName = getFileDirectory() + "level" + std::to_string(iLevel) + "_";
      std::string fileExtension = ".bin";
      cfgUtil::readBinary<int>(fileRootName + "matAssignments" + fileExtension, m_materialAssignments[iLevel]);

      /*
      Vector3f forceAxis;
      if (m_readFullDeformation)
      {
        std::string baseFileName = "StressDeformation_" + std::to_string(level) + "_" + std::to_string(m_blockSize);
        std::string stressStrainFileName = fileDirectory + baseFileName + "_x.txt";
        std::vector<std::vector<float> > stresses;
        std::vector<std::vector<std::vector<float> > > x;
        int n[3];
        ResOk = readData(stressStrainFileName,  forceAxis, m_materialAssignments[iLevel], stresses, x, n);
      }
      else
      {
        std::string baseFileName = "StressStrain_" + std::to_string(level) + "_" + std::to_string(m_blockSize);
        std::string stressStrainFileName = fileDirectory + baseFileName + "_x.txt";
        std::vector<std::vector<float> > stresses, strains;
        ResOk = readData(stressStrainFileName,  forceAxis, m_materialAssignments[iLevel], stresses, strains);
      }*/ 
    }
    materialAssignment = m_materialAssignments[iLevel][indComb];
  }
  else
  {
    std::string stressStrainFilesDir[2];
    stressStrainFilesDir[0] =  fileDirectory + "x_level" + std::to_string(level) + "//";
    stressStrainFilesDir[1] =  fileDirectory + "y_level" + std::to_string(level) + "//";
    std::string stressStrainBaseFileName = "StressStrain_" + std::to_string(level) + "_" + std::to_string(m_blockSize);

    std::string sampleFile = stressStrainFilesDir[0] + stressStrainBaseFileName + "_" + std::to_string(indComb) + ".txt"; 
    Vector3f forceAxis;
    
    std::vector<cfgScalar> stresses, strains;
    ResOk = readData(sampleFile, forceAxis, materialAssignment, stresses, strains);
  }
#endif 

  int level = m_levels[iLevel];
  int N[3] = {level,level,level};

  const std::vector<int> & materialAssignment = m_materialAssignments[iLevel][iCombIndex];

  bool isManifold = cfgMaterialUtilities::isStructureManifold(N[0], N[1], N[2], materialAssignment, N[0], N[1], N[2]);

  int blockSize = 3;
  int c = blockSize;
  int n[3] = {c*N[0], c*N[1], c*N[2]};

  std::vector<int> cellMaterials;
  cellMaterials = materialAssignment;
  if (m_dim==2)
  {
    repMaterialAssignment(N[0], N[1], materialAssignment, blockSize, blockSize, 1, cellMaterials);
    n[2] = 1;
  }
  else
  {
    repMaterialAssignment(N[0], N[1], N[2], materialAssignment, blockSize, blockSize, blockSize, 1, cellMaterials);
  }
  ElementRegGrid * em = new ElementRegGrid(n[0], n[1], n[2]);
  em->me = cellMaterials;

  QSharedPointer<ElementMesh> elementMesh(em);
  return elementMesh;
}

void exProject::setLevelVisibility(int ilevel, bool isVisible) 
{
  int visibility = (isVisible? 1: 0);
  bool updateView = (m_levelVisibility[ilevel] != visibility);
  m_levelVisibility[ilevel] = isVisible;
  if (updateView)
  {
    emit levelVisibilityModified();
  }
}

bool exProject::getLevelVisibility(int ilevel)
{
  assert(ilevel>=0 && ilevel<m_levelVisibility.size());
  return m_levelVisibility[ilevel];
}







