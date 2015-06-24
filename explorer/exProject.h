#ifndef exProject_h
#define exProject_h

#include "cfgDefs.h"

#include <QObject>
#include <QSharedPointer>

class ElementMesh;

class exProject: public QObject
{
  Q_OBJECT

public:
  exProject(int iDim);
  virtual ~exProject();

  void setFileDirectory(const std::string &iFileDirectory) {m_FileDirectory = iFileDirectory;}
  void setFileSubdirectories(const std::string &iFileDirectoryX, const std::string &iFileDirectoryY);

  bool getMaterialParameters(int iLevel, std::vector<cfgScalar> &oPhysicalParameters, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<std::vector<int> > &oMaterialAssignment);

  QSharedPointer<ElementMesh> computeElementMesh(int iCombIndex, int iLevel);

  void setPickedStructure(int iStructureIndex, int iStructureLevel, const std::vector<int> &iBaseMaterials);
  int getPickedStructureIndex() const {return m_pickedStructureIndex;}
  int getPickedStructureLevel() const {return m_pickedStructureLevel;}
  const std::vector<int> & getPickedStructureBaseMaterials() const {return m_pickedStructureBaseMaterials;}

signals:
  void pickedStructureModified();

private:
  bool readBaseMaterialStructures(int iLevel, std::vector<std::vector<int> > &oBaseMaterialStructures);

private:
  int m_blockSize;
  int m_dim;

  int m_pickedStructureIndex;
  int m_pickedStructureLevel;
  std::vector<int> m_pickedStructureBaseMaterials;

  std::string m_FileDirectory;
  std::string m_FileSubdirectories;
};

#endif

