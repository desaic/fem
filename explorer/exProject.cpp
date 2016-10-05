#include "exProject.h"

#include "ElementRegGrid.hpp"
#include "cfgDefs.h"
#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;
#include <cfgUtilities.h>
using namespace cfgUtil;

#include <assert.h>

#include <QStringList>

#include "ExplorerUtilities.h"
#include "FamilyExtractor.h"
#include "MicrostructureSet.h"
#include "FamilyGenerator.h"

exProject::exProject()
{
  m_blockSize = 1;
  m_dim = 3;
  m_readSingleFile = true;
  m_readFullDeformation = true;

  m_paramsToVisualize.push_back(EXType);
  m_paramsToVisualize.push_back(NuXYType);
  m_paramsToVisualize.push_back(DensityType);
  m_paramsToVisualize.push_back(MuXYType);

  m_activeTool = ex::None;

  m_selectedGamutIndex = -1;
}

exProject::~exProject()
{
}

void exProject::setDimension(int iDim)
{
  Q_ASSERT(iDim==2 || iDim==3);
  m_dim = iDim;
}

void exProject::setParameterToVisualize(int indParam, MaterialParameterType iType)
{
  assert(indParam>=0 && indParam < 4);
  m_paramsToVisualize[indParam] = iType;
}

void exProject::updateParameterToVisualize(int indParam, MaterialParameterType iType)
{
  setParameterToVisualize(indParam, iType);
  emit paramsToVisualizeModified();
}

exProject::MaterialParameterType exProject::getParameterToVisualize(int indParam)
{
  assert(indParam>=0 && indParam < 4);
  return m_paramsToVisualize[indParam];
}

void exProject::setPickedReducedPoint(int iPointIndex) 
{
  std::cout << "exProject::setPickedReducedPoint, iPointIndex = " << iPointIndex << " npoint =  "<< m_microstructuresIndices.size() << std::endl;
  assert(iPointIndex>=0 && iPointIndex<m_microstructuresIndices.size());
  if (iPointIndex>=0 && iPointIndex<m_microstructuresIndices.size())
  {
    m_pickedStructureIndex = m_microstructuresIndices[iPointIndex];
    m_pickedStructureLevel = 0;

    emit pickedStructureModified();
    emit pickedReducedPointModified();
  }
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

    std::vector<cfgScalar> physicalParameters, elasticityTensors, thermalExpansionCoeffs, ultimateStrengths;
    std::vector<std::vector<int> > baseMaterials, materialAssignments;
    std::vector<std::vector<cfgScalar> > materialDistributions;
    resOk = cfgUtil::readBinary<float>(fileRootName + "params" + fileExtension, physicalParameters);
    //if (resOk)
    //  resOk = cfgUtil::readBinary<int>(fileRootName + "baseMat" + fileExtension, baseMaterials);
    if (resOk)
    {
      cfgUtil::readBinary<int>(fileRootName + "matAssignments" + fileExtension, materialAssignments);
      cfgUtil::readBinary<cfgScalar>(fileRootName + "matDistributions" + fileExtension, materialDistributions);
      cfgUtil::readBinary<cfgScalar>(fileRootName + "elasticityTensors" + fileExtension, elasticityTensors);
      cfgUtil::readBinary<cfgScalar>(fileRootName + "thermalExpansionCoefficients" + fileExtension, thermalExpansionCoeffs);
      cfgUtil::readBinary<cfgScalar>(fileRootName + "ultimateStrengths" + fileExtension, ultimateStrengths);
    }

    if (resOk)
    {
      oLevel = level;

      QString str = iFileName.section("level", 0, 0);
      oLabel = stringRoot;
      oLabel.chop(1);
      oLabel.remove(str);

      addLevelData(level, baseMaterials, materialAssignments, materialDistributions, physicalParameters, elasticityTensors, thermalExpansionCoeffs, ultimateStrengths);
    }
  }
  std::cout << (resOk? "ok": "ko") << std::endl;
  return resOk;
}

void exProject::addLevelData(int ilevel, const std::vector<std::vector<int> > &iBaseMaterials, const std::vector<std::vector<int> > &iMaterialAssignments, const std::vector<std::vector<cfgScalar> > &iMaterialDistributions, 
                             const std::vector<cfgScalar> &iPhysicalParameters, const std::vector<cfgScalar> &iElasticityTensors, const std::vector<cfgScalar> &iThermalExpansionCoeffs, const std::vector<cfgScalar> &iUltimateStrengths)
{
  m_levels.push_back(ilevel);
  m_baseMaterials.push_back(iBaseMaterials);
  m_materialAssignments.push_back(iMaterialAssignments);
  m_materialDistributions.push_back(iMaterialDistributions);
  m_physicalParametersPerLevel.push_back(iPhysicalParameters);
  m_elasticityTensors.push_back(iElasticityTensors);
  m_thermalExpansionCoeffs.push_back(iThermalExpansionCoeffs);
  m_ultimateStrengths.push_back(iUltimateStrengths);

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

  int nlevel = (int)m_levels.size()+1;

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
    Vector3S forceAxis;
    
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

  std::vector<int> materialAssignment;
  if (m_materialAssignments[iLevel].size()>0)
  {
    materialAssignment = m_materialAssignments[iLevel][iCombIndex];
  }
  else if (m_materialDistributions[iLevel].size()>0)
  {
    materialAssignment = convertVec<float, int>(mult(m_materialDistributions[iLevel][iCombIndex], 255.f));
  }

  //bool isManifold = cfgMaterialUtilities::isStructureManifold(N[0], N[1], N[2], materialAssignment, true, 1, 1, 1);
  //std::cout << "isManifold = " << (isManifold? "true": "false") << std::endl;

  int blockSize = 3;
  if (m_dim==3)
  {
    blockSize = 1;
  }
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

void exProject::runFamilyGenerator()
{
  int n[3] = {16, 16, m_dim==2?1: 16};
  FamilyGenerator * familyGenerator = FamilyGenerator::createOperator();
  familyGenerator->setMicrostructureSize(m_dim, n);
  familyGenerator->run();
  const std::vector<std::vector<int> > & matAssignments = familyGenerator->getMicrostructures();
  std::vector<cfgScalar> parameters;
  //const std::vector<cfgScalar> & parameters = familyGenerator->getParameters();
  familyGenerator->getSize(n);
  ExplorerUtilities::writeMicrostructures("..//..//Output//", n[0], matAssignments, parameters, "family");

  SAFE_DELETE(familyGenerator);
}

void exProject::runFamilyExtractor(int iOption)
{
  if (iOption==1)
  {
    FamilyExtractor * familyExtractor = FamilyExtractor::createOperator();
    familyExtractor->setOption(iOption);
    familyExtractor->run();
    m_microstructuresReducedCoords = familyExtractor->getReducedCoordinates();
    m_microstructuresIndices = familyExtractor->getMicrostructureIndices();
    m_selectedGamutIndex = 0;
    SAFE_DELETE(familyExtractor);
  }
  else
  {
    int ind = m_selectedGamutIndex;
    if (ind<0)
    {
      m_selectedGamutIndex = 0;
      std::string fileNameIndices = "..//..//Output//microstructureIndices.txt";
      cfgUtil::readVectorFromFile<int>(fileNameIndices, m_selectedPointIndices); 
      ind = m_selectedGamutIndex;
    }
    if (ind>=0)
    {
      FamilyExtractor * familyExtractor = FamilyExtractor::createOperator();

      int microstructureSize[3] = {m_levels[ind], m_levels[ind], m_levels[ind]};
      int paramDim = (int)(m_physicalParametersPerLevel[ind].size()/m_materialAssignments[ind].size());

      MicrostructureSet microstructures;
      microstructures.setMicrostructures(&m_materialAssignments[ind], microstructureSize, &m_physicalParametersPerLevel[ind], paramDim);

      familyExtractor->setMicrostructures(&microstructures, &m_selectedPointIndices);

      familyExtractor->setOption(iOption);

      familyExtractor->run();

      m_microstructuresReducedCoords = familyExtractor->getReducedCoordinates();
      m_microstructuresIndices = familyExtractor->getMicrostructureIndices();

      SAFE_DELETE(familyExtractor);
    }
  }
}








