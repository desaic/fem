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
  exProject();
  virtual ~exProject();

  void setFileDirectory(const std::string &iFileDirectory) {m_FileDirectory = iFileDirectory;}
  void setFileSubdirectories(const std::string &iFileDirectoryX, const std::string &iFileDirectoryY);

  bool getMaterialParameters(std::vector<cfgScalar> &oPhysicalParameters, int iLevel);

  QSharedPointer<ElementMesh> computeElementMesh(int iCombIndex, int iLevel);

  void setPickedStructure(int iStructureIndex, int iStructureLevel);
  int getPickedStructureIndex() const {return m_pickedStructureIndex;}
  int getPickedStructureLevel() const {return m_pickedStructureLevel;}

signals:
  void pickedStructureModified();

private:
  bool readBaseMaterialStructures(int iLevel, std::vector<std::vector<int> > &oBaseMaterialStructures);

private:
  int m_blockSize;

  int m_pickedStructureIndex;
  int m_pickedStructureLevel;

  std::string m_FileDirectory;
  std::string m_FileSubdirectories;
};

#endif

