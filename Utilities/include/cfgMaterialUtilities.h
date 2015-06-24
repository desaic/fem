#ifndef cfgMaterialUtilities_h
#define cfgMaterialUtilities_h

#include <string>
#include <map>

#include "cfgDefs.h"

#include <Vector3f.h>

class ElementRegGrid;
class ElementRegGrid2D;

namespace cfgMaterialUtilities
{
  // I/O
  // ---
  bool readMaterialCombinations(const std::string iFileName, std::vector<std::vector<int> > &oMaterials);
  bool writeMaterialCombinations(const std::string iFileName, const std::vector<std::vector<int> > &iMaterials);
  bool readData(const std::string &iFileName,  Vector3f &oForceAxis, std::vector<int> &oMaterials, std::vector<float> &oStresses, std::vector<float> &oStrains);
  bool writeData(const std::string &iFileName, const Vector3f &iForceAxis, const std::vector<int> &iMaterials, const std::vector<float> &iStresses, const std::vector<float> &iStrains);
  bool computeMaterialParameters(int idim, const std::string &iMaterialFile, const std::string iStressStrainFilesDirectories[2], const std::string iStressStrainBaseFileName, int iNbFiles, 
                                 std::vector<cfgScalar> &oPhysicalParameters, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<std::vector<int> > &oMaterialAssignments,
                                 std::vector<std::vector<int> > &oMaterialAssignmentsOneCell);

  cfgScalar computeYoungModulus(const std::vector<cfgScalar> &iStrains, const std::vector<cfgScalar>  &iStresses);

  // Hexahedral mesh utilities
  // -------------------------
  int getElementIndex(int nx, int ny, int nz, int ii , int jj, int kk);
  void getMaterialAssignment(int nx, int ny, int nz, const std::vector<int> &iMaterialCombIndices, int nX, int nY, int nZ, const std::vector<std::vector<int> > &iBaseCellMaterials, std::vector<int> &oMaterials);

  //0: Left, 1: Right, 2: Bottom, 3:Top, 4: Back, 5: Front
  void getSideVoxels(int iSide, int nx, int ny, int nz, std::vector<int> &oSideVoxelIndices);
  void getSideVertices(int iSide, const ElementRegGrid * iElementGrid, std::vector<int> &oElemIndices, std::vector<std::vector<int> > &oFaceVertexIndices);
  void getExternalVertices(ElementRegGrid * iElementGrid, std::vector<int> &oElemIndices, std::vector<std::vector<int> > &oFaceVertexIndices, std::vector<Vector3f> &oNormals);
  void getVertexValences(const ElementRegGrid * iElementGrid, const std::vector<int> &iElemIndices, const std::vector<std::vector<int> > &iFvIndices, std::map<int,int> &oind2Valence);

  // Hexahedral mesh utilities 2D
  // ----------------------------
  void getSideVertices(int iSide, const ElementRegGrid2D * iElementGrid, std::vector<int> &oElemIndices, std::vector<std::vector<int> > &oFaceVertexIndices);
  void getVertexValences(const ElementRegGrid2D * iElementGrid, const std::vector<int> &iElemIndices, const std::vector<std::vector<int> > &iEvIndices, std::map<int,int> &oind2Valence);
  void getMaterialAssignment(int nx, int ny, const std::vector<int> &iMaterialCombIndices, int nX, int nY, const std::vector<std::vector<int> > &iBaseCellMaterials, std::vector<int> &oMaterials);

  // Point cloud utitlities
  // ----------------------
  void computeConvexHull(const std::vector<float> &iPoints, int iDim, std::vector<int> &oConvexHullVertices);
  void computeDelaundayTriangulation(const std::vector<float> &iPoints, int iDim, std::vector<int> &oFaces, std::vector<int> &oBoundaryVertices, std::vector<int> &oBoundaryFaces);

  // Triangle mesh utilities
  //------------------------
  float getMeanEdgeLength(const std::vector<float> &iX,  const std::vector<int> &iIndexArray, int iDim);
  Vector3f getMeanEdgeLengthPerAxis(const std::vector<float> &iX,  const std::vector<int> &iIndexArray, int iDim);
  float getMinEdgeLength(const std::vector<float> &iX,  const std::vector<int> &iIndexArray, int iDim);
  void getEdgesFromTriFaceIndexArray(const std::vector<int> &iTriIndexArray, std::vector<int> &oEdgeIndexArray);

  // Tet mesh utilities
  // ------------------
  void getEdgesFromTetFaceIndexArray(const std::vector<int> &iTetIndexArray, std::vector<int> &oEdgeIndexArray);

  // Conversion between types
  // -----------------------
  Vector3f getVector3f(int indVertex, const std::vector<float> &iPoints);
  Vector2f getVector2f(int indVertex, const std::vector<float> &iPoints);
};

#endif 

