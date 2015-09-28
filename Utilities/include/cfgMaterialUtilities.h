#ifndef cfgMaterialUtilities_h
#define cfgMaterialUtilities_h

#include <string>
#include <map>

#include "cfgDefs.h"

#include <Vector3f.h>

namespace cfgMaterialUtilities
{
  // I/O
  // ---
  bool readMaterialCombinations(const std::string iFileName, std::vector<std::vector<int> > &oMaterials);
  bool writeMaterialCombinations(const std::string iFileName, const std::vector<std::vector<int> > &iMaterials);
  bool readData(const std::string &iFileName,  Vector3f &oForceAxis, std::vector<int> &oMaterials, std::vector<float> &oStresses, std::vector<float> &oStrains);
  bool writeData(const std::string &iFileName, const Vector3f &iForceAxis, const std::vector<int> &iMaterials, const std::vector<float> &iStresses, const std::vector<float> &iStrains); 
  bool writeData(const std::string &iFileName, const Vector3f &iForceAxis, const std::vector<std::vector<int> > &iMaterials, const std::vector<std::vector<float> > &iStresses, const std::vector<std::vector<float> > &iStrains);
  bool writeData(const std::string &iFileName, const Vector3f &iForceAxis, const std::vector<std::vector<int> > &iMaterials, const std::vector<std::vector<float> > &iStresses, const std::vector<std::vector<std::vector<float> > > &ix, int iN[3]);
  bool writeDataBinary(const std::string &iFileName, const Vector3f &iForceAxis, const std::vector<std::vector<int> > &iMaterials, const std::vector<std::vector<float> > &iStresses, const std::vector<std::vector<std::vector<float> > > &ix, int iN[3]);
  bool concatenateData(const std::string &iInputBaseFileName, int iNbFiles, std::string &iOutputFileName);
  bool readData(const std::string &iFileName,  Vector3f &oForceAxis, std::vector<std::vector<int> > &oMaterials, std::vector<std::vector<float> > &oStresses, std::vector<std::vector<float> > &oStrains);
  bool readData(const std::string &iFileName,  Vector3f &oForceAxis, std::vector<std::vector<int> > &oMaterials, std::vector<std::vector<float> > &oStresses, std::vector<std::vector<std::vector<float> > > &ox, int oN[3]);
  bool readDataBinary(const std::string &iFileName,  Vector3f &oForceAxis, std::vector<std::vector<int> > &oMaterials, std::vector<std::vector<float> > &oStresses, std::vector<std::vector<std::vector<float> > > &ox, int oN[3]);
  
  bool computeMaterialParameters(const std::string &iMaterialFile, const std::string iStressStrainFilesDirectories[2], const std::string iStressStrainBaseFileName, int iNbFiles, 
                                 std::vector<cfgScalar> &oPhysicalParameters, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<std::vector<int> > &oMaterialAssignmentsOneCell);
  bool computeMaterialParameters(const std::string &iMaterialFile, const std::string iStressStrainFileNames[2],
                                 std::vector<cfgScalar> &oPhysicalParameters, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<std::vector<int> > &oMaterialAssignmentsOneCell);
  bool computeMaterialParametersFromDeformations(const std::string &iMaterialFile, const std::string iStressDeformationFileNames[2], bool iBinary,
                                                 std::vector<cfgScalar> &oPhysicalParameters, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<std::vector<int> > &oMaterialAssignmentsOneCell);

   bool computeMaterialParameters(const std::vector<std::vector<int> > &iMaterials, const std::vector<std::vector<int> > &iBaseMaterialStructures,
                                 const std::vector<std::vector<cfgScalar> > iStresses[2], const std::vector<std::vector<cfgScalar> > iStrains[2], std::vector<cfgScalar> &oPhysicalParameters);

  bool computeMaterialParameters(const std::vector<std::vector<int> > &iMaterials, const std::vector<std::vector<int> > &iBaseMaterialStructures,
                                 const std::vector<std::vector<cfgScalar> > iStresses[2], const std::vector<std::vector<cfgScalar> > iStrains[2][2], std::vector<cfgScalar> &oPhysicalParameters);

  bool getMaterialAssignment(int iDim, std::string iMaterialFile, const std::vector<int> &iMaterialsCombination, int iBlockRep, int iNbSubdiv, std::vector<int> &oCellMaterials);

  cfgScalar computeYoungModulus(const std::vector<cfgScalar> &iStrains, const std::vector<cfgScalar>  &iStresses);
  cfgScalar computeMaterialRatio(const std::vector<int> &iMaterialAssignment, const std::vector<std::vector<int> > &iBaseMaterialStructures);
  cfgScalar computePoissonRatio(const std::vector<cfgScalar> &iStrainsAxialDir, const std::vector<cfgScalar> &iStrainsTransversalDir);
  cfgScalar computeVonMisesStress(std::vector<cfgScalar> &iStresses);  // iStresses = [sigma_xx sigma_yy sigma_zz sigma_xy sigma_yz sigma_zx]
  cfgScalar computeVonMisesStress2D(std::vector<cfgScalar> &iStresses); //iStresses = [sigma_xx sigma_yy sigma_xy]

  // Hexahedral mesh utilities
  // -------------------------
  void getMaterialAssignment(int nx, int ny, int nz, const std::vector<int> &iMaterialCombIndices, int nX, int nY, int nZ, const std::vector<std::vector<int> > &iBaseCellMaterials, int repX, int repY, int repZ, int iNbSubdivisions, 
                             std::vector<int> &oMaterials);

  // Hexahedral mesh utilities 2D
  // ----------------------------
  int getGridToVectorIndex(int i, int j, int nx, int ny);
  void getVectorIndexToGrid(int iIndex, int nx, int ny, int &oIndex_i, int &oIndex_j);

  //0: Left, 1: Right, 2: Bottom, 3:Top
  void getSideElements(int iSide, int nx, int ny, std::vector<int> &oElementIndices);
  void getSideVertices(int iSide, int nx, int ny, std::vector<int> &oVertexIndices);
  void getBoundaryCellIndices(int nx, int ny, std::vector<int> &oCellIndices);
  void getMaterialAssignment(int nx, int ny, const std::vector<int> &iMaterialCombIndices, int nX, int nY, const std::vector<std::vector<int> > &iBaseCellMaterials, int repX, int repY, int iNbSubdivisions, std::vector<int> &oMaterials);

  void repMaterialAssignment(int nx, int ny, const std::vector<int> &iMaterialCombIndices, int iNbSubdivions, int repX, int repY, std::vector<int> &oMaterials);
  bool isStructureManifold(int nx, int ny, const std::vector<int> &iMaterials, int nX, int nY, int iNbCellsToCheck);
  bool isStructureManifold(int nx, int ny, const std::vector<int> &iMaterials);
  void updateMaterialSubstructure(std::vector<int> &ioMaterialAsssignment, int nx, int ny, const std::vector<int> &iSubMaterialStructure, int Nx, int Ny, int iSubMatStructLocation);
  
  void getLayer(int Nx, int NY, const std::vector<int> &iStructureElements, int iIndex, int iDim, std::vector<int> &oNewLayerElements);
  void insertLayer(int Nx, int Ny, const std::vector<int> &iStructureElements, std::vector<int> &iNewLayerElements, int iIndex, int iDim, std::vector<int> &oNewStructureElements);

  cfgScalar computeStrain(const std::vector<cfgScalar> &ix, int nx, int ny, int iAxis);
 
  int getMaterialSignature(int nx, int ny, const std::vector<int> &icellMaterials, const std::vector<int> &iIndicesToUse);
  void clusterByBoundary(int nx, int ny, const std::vector<std::vector<int> > &icellMaterials, int iNbMat, std::vector<std::vector<int> > &oClusters);

  // Point cloud utitlities
  // ----------------------
  void computeConvexHull(const std::vector<float> &iPoints, int iDim, std::vector<int> &oConvexHullVertices);
  void computeDelaundayTriangulation(const std::vector<float> &iPoints, int iDim, std::vector<int> &oFaces, std::vector<int> &oBoundaryVertices, std::vector<int> &oBoundaryFaces, std::vector<float> *oDistancesToBoundary=NULL);
  void getClosestPoints(const std::vector<float> &iPoints, int iDim, std::vector<int> &iRefPointIndices, float iRange, std::vector<int> &oPoints);
  void getFurthestPointsGreedy(int iOutputNbPoints, const std::vector<cfgScalar> &iPoints, int iDim, std::vector<int> &oPointIndices);
  void getKMeans(int iNbIterations, int iNbClusters, const std::vector<cfgScalar> &iPoints, int iDim, std::vector<std::vector<int> > &oClusters, std::vector<int> *oCenters=NULL);

  void getBoundingBox(const std::vector<cfgScalar> &iPoints, int iDim, std::vector<cfgScalar> oBox[2]);
  void rescaleData(std::vector<cfgScalar> &ioPoints, int iDim,  const std::vector<cfgScalar> &iTargetBoxLengths);

  // Triangle mesh utilities
  //------------------------
  float getMeanEdgeLength(const std::vector<float> &iX,  const std::vector<int> &iIndexArray, int iDim);
  Vector3f getMeanEdgeLengthPerAxis(const std::vector<float> &iX,  const std::vector<int> &iIndexArray, int iDim);
  Vector3f getMedianEdgeLengthPerAxis(const std::vector<float> &iX,  const std::vector<int> &iIndexArray, int iDim);
  float getMinEdgeLength(const std::vector<float> &iX,  const std::vector<int> &iIndexArray, int iDim);
  void getEdgesFromTriFaceIndexArray(const std::vector<int> &iTriIndexArray, std::vector<int> &oEdgeIndexArray);
  void sampleMesh(const std::vector<float> &iPoints, const std::vector<int> &iTriIndexArray, int iNbSamplesPerFace, std::vector<float> &oPoints);
  void computeDistancesToPointCloud(const std::vector<float> &iPoints, const std::vector<float> &iPointCloud, std::vector<float> &oDistances);

  // Tet mesh utilities
  // ------------------
  void getEdgesFromTetFaceIndexArray(const std::vector<int> &iTetIndexArray, std::vector<int> &oEdgeIndexArray);

  // Conversion between types
  // -----------------------
  Vector3f getVector3f(int indVertex, const std::vector<float> &iPoints);
  Vector2f getVector2f(int indVertex, const std::vector<float> &iPoints);
  Vector2S getVector2S(int indVertex, const std::vector<cfgScalar> &iPoints);

  std::vector<float> toVectorFloat(const std::vector<Vector3f> &iPoints);
  std::vector<float> toVectorFloat(const std::vector<Vector2f> &iPoints);

  std::vector<cfgScalar> toVectorScalar(const std::vector<Vector3S> &iPoints);
  std::vector<cfgScalar> toVectorScalar(const std::vector<Vector2S> &iPoints);

  std::vector<Vector2f> toVector2f(const std::vector<float> &iPoints);
  std::vector<Vector2S> toVector2S(const std::vector<cfgScalar> &iPoints);
  std::vector<Vector3f> toVector3f(const std::vector<float> &iPoints);
};

#endif 

