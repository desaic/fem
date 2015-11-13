#ifndef exProject_h
#define exProject_h

#include "cfgDefs.h"

#include <QObject>
#include <QSharedPointer>

class ElementMesh;

#define MAX_LEVELS 10000

class exProject: public QObject
{
  Q_OBJECT

public:
  exProject(int iDim);
  virtual ~exProject();

  bool loadFile(const QString &iFileName, int &oLevel, QString &oLabel);

  void setFileDirectory(const std::string &iFileDirectory) {m_FileDirectory = iFileDirectory;}
  void setFileSubdirectories(const std::string &iFileDirectoryX, const std::string &iFileDirectoryY);
  const std::string & getFileDirectory() {return m_FileDirectory;}

  void setMaterialAssignment(int iLevel, const std::vector<std::vector<int> > &iMaterialAssignment);
  bool getMaterialParameters(int iLevel, std::vector<cfgScalar> &oPhysicalParameters, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<std::vector<int> > &oMaterialAssignmentPerCell);

  QSharedPointer<ElementMesh> computeElementMesh(int iCombIndex, int iLevel);
  QSharedPointer<ElementMesh> computeElementMeshIncr(int iCombIndex, int iLevel);

  void setPickedStructure(int iStructureIndex, int iStructureLevel);
  int getPickedStructureIndex() const {return m_pickedStructureIndex;}
  int getPickedStructureLevel() const {return m_pickedStructureLevel;}

  void setLevelVisibility(int ilevel, bool isVisible);
  bool getLevelVisibility(int ilevel);

  const std::vector<std::vector<cfgScalar> > & getPhysicalParameters() {return m_physicalParametersPerLevel;}
  const std::vector<std::vector<std::vector<int> > > & getMaterialAssignments() {return m_materialAssignments;}
  const std::vector<int>  & getLevels() {return m_levels;}

signals:
  void pickedStructureModified();
  void levelVisibilityModified();
  void levelsModified();

private:
  bool readBaseMaterialStructures(int iLevel, std::vector<std::vector<int> > &oBaseMaterialStructures);

  void addLevelData(int ilevel, const std::vector<std::vector<int> > &iBaseMaterials, const std::vector<std::vector<int> > &iMaterialAssignments, std::vector<cfgScalar> &iPhysicalParameters);

private:
  int m_blockSize;
  int m_dim;

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
  std::vector<std::vector<cfgScalar> > m_physicalParametersPerLevel;

  std::vector<int> m_levels;
};

#endif

