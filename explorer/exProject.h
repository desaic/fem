#ifndef exProject_h
#define exProject_h

#include "cfgDefs.h"

#include <QObject>
#include <QSharedPointer>

class ElementMesh;

#define MAX_LEVELS 40

class exProject: public QObject
{
  Q_OBJECT

public:
  exProject(int iDim);
  virtual ~exProject();

  void setFileDirectory(const std::string &iFileDirectory) {m_FileDirectory = iFileDirectory;}
  void setFileSubdirectories(const std::string &iFileDirectoryX, const std::string &iFileDirectoryY);
  const std::string & getFileDirectory() {return m_FileDirectory;}

  void setMaterialAssignment(int iLevel, const std::vector<std::vector<int> > &iMaterialAssignment);
  bool getMaterialParameters(int iLevel, std::vector<cfgScalar> &oPhysicalParameters, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<std::vector<int> > &oMaterialAssignmentPerCell);

  QSharedPointer<ElementMesh> computeElementMesh(int iCombIndex, int iLevel);
  QSharedPointer<ElementMesh> computeElementMeshIncr(int iCombIndex, int iLevel);

  void setPickedStructure(int iStructureIndex, int iStructureLevel, const std::vector<int> &iBaseMaterials);
  int getPickedStructureIndex() const {return m_pickedStructureIndex;}
  int getPickedStructureLevel() const {return m_pickedStructureLevel;}
  const std::vector<int> & getPickedStructureBaseMaterials() const {return m_pickedStructureBaseMaterials;}

  int getMaxLevel() {return m_maxLevel;}

  void setLevelVisibility(int ilevel, bool isVisible);
  bool getLevelVisibility(int ilevel);

signals:
  void pickedStructureModified();
  void levelVisibilityModified();

private:
  bool readBaseMaterialStructures(int iLevel, std::vector<std::vector<int> > &oBaseMaterialStructures);

private:
  int m_blockSize;
  int m_dim;
  int m_maxLevel;

  int m_pickedStructureIndex;
  int m_pickedStructureLevel;
  std::vector<int> m_pickedStructureBaseMaterials;

  bool m_readSingleFile;
  bool m_readFullDeformation;
  std::vector<std::vector<int> > m_materialAssignments[MAX_LEVELS];
  std::string m_FileDirectory;
  std::string m_FileSubdirectories;

  std::vector<int> m_levelVisibility;
};

#endif

