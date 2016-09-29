#ifndef exProject_h
#define exProject_h

#include "cfgDefs.h"

#include <QObject>
#include <QSharedPointer>
#include "ExplorerDefs.h"

class ElementMesh;

#define MAX_LEVELS 10000

class exProject: public QObject
{
  Q_OBJECT

public:
  enum MicrostructureType
  {
    CubicType=0,
    OrthotropicType
  };

  enum MaterialParameterType
  {
    EXType=0,
    EYType,
    EZType,
    NuXYType,
    NuXZType,
    NuYZType,
    MuXYType,
    MuXZType,
    MuYZType,
    DensityType,
    StrengthType,
    UndefinedType,
  };

public:
  exProject();
  virtual ~exProject();

  void setDimension(int iDim);
  void setType(MicrostructureType iType) {m_type = iType;}
  void setParameterToVisualize(int indParam, MaterialParameterType iType);
  void updateParameterToVisualize(int indParam, MaterialParameterType iType);

  bool loadFile(const QString &iFileName, int &oLevel, QString &oLabel);

  void setFileDirectory(const std::string &iFileDirectory) {m_FileDirectory = iFileDirectory;}
  void setFileSubdirectories(const std::string &iFileDirectoryX, const std::string &iFileDirectoryY);
  const std::string & getFileDirectory() {return m_FileDirectory;}

  void setMaterialAssignment(int iLevel, const std::vector<std::vector<int> > &iMaterialAssignment);
  bool getMaterialParameters(int iLevel, std::vector<cfgScalar> &oPhysicalParameters, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<std::vector<int> > &oMaterialAssignmentPerCell);

  QSharedPointer<ElementMesh> computeElementMesh(int iCombIndex, int iLevel);
  QSharedPointer<ElementMesh> computeElementMeshIncr(int iCombIndex, int iLevel);

  void setPickedReducedPoint(int iPointIndex);

  void clearSelectedPoints() {m_selectedPointIndices.clear(); m_selectedGamutIndex=-1;}
  void setSelectedPoints(const std::vector<int> &iPointIndices, int iGamutIndex) {m_selectedPointIndices = iPointIndices; m_selectedGamutIndex = iGamutIndex;}
  int getSelectedGamutIndex() {return m_selectedGamutIndex;}

  void setPickedStructure(int iStructureIndex, int iStructureLevel);
  int getPickedStructureIndex() const {return m_pickedStructureIndex;}
  int getPickedStructureLevel() const {return m_pickedStructureLevel;}

  void setLevelVisibility(int ilevel, bool isVisible);
  bool getLevelVisibility(int ilevel);

  void setDistancesToBoundary(const std::vector<std::vector<float> > &iDistances) { m_distancesToBoundary = iDistances;};
  const std::vector<std::vector<float> > & getDistancesToBoundary() {return m_distancesToBoundary;}

  void setPhysicalParameters(const std::vector<std::vector<cfgScalar> > &iParams) {m_physicalParametersPerLevel = iParams;}
  const std::vector<std::vector<cfgScalar> > & getPhysicalParameters() {return m_physicalParametersPerLevel;}
  const std::vector<std::vector<cfgScalar> > & getElasticityTensors() {return m_elasticityTensors;}
  const std::vector<std::vector<cfgScalar> > & getThermalExpansionCoeffs() {return m_thermalExpansionCoeffs;}
  const std::vector<std::vector<cfgScalar> > & getUltimateStrengths() {return m_ultimateStrengths;}
  const std::vector<std::vector<std::vector<int> > > & getMaterialAssignments() {return m_materialAssignments;}
  const std::vector<int>  & getLevels() {return m_levels;}

  int getDim() {return m_dim;}
  MicrostructureType getType() {return m_type;}
  MaterialParameterType getParameterToVisualize(int indParam);

  void runFamilyGenerator();
  void runFamilyExtractor(int iOption);
  std::vector<cfgScalar> & getMicrostructuresReducedCoordinates() {return m_microstructuresReducedCoords;}
  std::vector<int> & getMicrostructuresIndices() {return m_microstructuresIndices;}

  void setActiveTool(ex::ToolType iType) {m_activeTool = iType;}
  ex::ToolType getActiveTool() {return m_activeTool;}

signals:
  void pickedStructureModified();
  void pickedReducedPointModified();
  void levelVisibilityModified();
  void levelsModified();
  void paramsToVisualizeModified();

private:
  bool readBaseMaterialStructures(int iLevel, std::vector<std::vector<int> > &oBaseMaterialStructures);

  void addLevelData(int ilevel, const std::vector<std::vector<int> > &iBaseMaterials, const std::vector<std::vector<int> > &iMaterialAssignments, const std::vector<std::vector<cfgScalar> > &iMaterialDistributions, 
                    const std::vector<cfgScalar> &iPhysicalParameters, const std::vector<cfgScalar> &iElasticityTensors, const std::vector<cfgScalar> &iThermalExpansionCoeffs, const std::vector<cfgScalar> &iUltimateStrengths);

private:
  int m_blockSize;
  int m_dim;
  MicrostructureType m_type;
  std::vector<MaterialParameterType> m_paramsToVisualize;

  int m_pickedStructureIndex;
  int m_pickedStructureLevel;

  bool m_readSingleFile;
  bool m_readFullDeformation;
  std::string m_FileDirectory;
  std::string m_FileSubdirectories;

  std::vector<int> m_levelVisibility;

  std::vector<std::vector<std::vector<int> > > m_baseMaterials;
  std::vector<std::map<std::vector<int>, int> > m_baseMaterials2Indices;
  std::vector<std::vector<std::vector<int> > > m_materialAssignments;
  std::vector<std::vector<std::vector<cfgScalar> > > m_materialDistributions;
  std::vector<std::vector<cfgScalar> > m_physicalParametersPerLevel;
  std::vector<std::vector<cfgScalar> > m_elasticityTensors;
  std::vector<std::vector<cfgScalar> > m_thermalExpansionCoeffs;
  std::vector<std::vector<cfgScalar> > m_ultimateStrengths;

  std::vector<int> m_levels;

  std::vector<std::vector<float> > m_distancesToBoundary;
  std::vector<int> m_microstructuresIndices;
  std::vector<cfgScalar> m_microstructuresReducedCoords;

  ex::ToolType m_activeTool;

  std::vector<int> m_selectedPointIndices;
  int m_selectedGamutIndex;
};

#endif

