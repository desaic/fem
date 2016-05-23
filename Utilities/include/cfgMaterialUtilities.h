#ifndef cfgMaterialUtilities_h
#define cfgMaterialUtilities_h
#include <Eigen/Dense>
#include <string>
#include <map>

#include "cfgDefs.h"

namespace cfgMaterialUtilities
{
  // I/O
  // ---
  bool readMaterialCombinations(const std::string iFileName, std::vector<std::vector<int> > &oMaterials);
  bool writeMaterialCombinations(const std::string iFileName, const std::vector<std::vector<int> > &iMaterials);
  bool readData(const std::string &iFileName, Eigen::Vector3f &oForceAxis, std::vector<int> &oMaterials, std::vector<float> &oStresses, std::vector<float> &oStrains);
  bool writeData(const std::string &iFileName, const Eigen::Vector3f &iForceAxis, const std::vector<int> &iMaterials, const std::vector<float> &iStresses, const std::vector<float> &iStrains);
  bool writeData(const std::string &iFileName, const Eigen::Vector3f &iForceAxis, const std::vector<std::vector<int> > &iMaterials, const std::vector<std::vector<float> > &iStresses, const std::vector<std::vector<float> > &iStrains);
  bool writeData(const std::string &iFileName, const Eigen::Vector3f &iForceAxis, const std::vector<std::vector<int> > &iMaterials, const std::vector<std::vector<float> > &iStresses, const std::vector<std::vector<std::vector<float> > > &ix, int iN[3]);
  bool writeDataBinary(const std::string &iFileName, const Eigen::Vector3f &iForceAxis, const std::vector<std::vector<int> > &iMaterials, const std::vector<std::vector<float> > &iStresses, const std::vector<std::vector<std::vector<float> > > &ix, int iN[3]);
  bool concatenateData(const std::string &iInputBaseFileName, int iNbFiles, std::string &iOutputFileName);
  bool readData(const std::string &iFileName, Eigen::Vector3f &oForceAxis, std::vector<std::vector<int> > &oMaterials, std::vector<std::vector<float> > &oStresses, std::vector<std::vector<float> > &oStrains);
  bool readData(const std::string &iFileName, Eigen::Vector3f &oForceAxis, std::vector<std::vector<int> > &oMaterials, std::vector<std::vector<float> > &oStresses, std::vector<std::vector<std::vector<float> > > &ox, int oN[3]);
  bool readDataBinary(const std::string &iFileName, Eigen::Vector3f &oForceAxis, std::vector<std::vector<int> > &oMaterials, std::vector<std::vector<float> > &oStresses, std::vector<std::vector<std::vector<float> > > &ox, int oN[3]);
  
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

   bool computeMaterialParameters(const std::vector<std::vector<int> > &iMaterials, const std::vector<std::vector<int> > &iBaseMaterialStructures,
                                 const std::vector<std::vector<cfgScalar> > iStresses[3], const std::vector<std::vector<cfgScalar> > iStrains[3][3], std::vector<cfgScalar> &oPhysicalParameters);

  bool getMaterialAssignment(int iDim, std::string iMaterialFile, const std::vector<int> &iMaterialsCombination, int iBlockRep, int iNbSubdiv, std::vector<int> &oCellMaterials);

  // Material properties
  // -------------------
  cfgScalar computeYoungModulus(cfgScalar iStrain, cfgScalar iStress);
  cfgScalar computeYoungModulus(const std::vector<cfgScalar> &iStrains, const std::vector<cfgScalar>  &iStresses);
  cfgScalar computeMaterialRatio(const std::vector<int> &iMaterialAssignment, const std::vector<std::vector<int> > &iBaseMaterialStructures);
  cfgScalar computeMaterialDensity(const std::vector<int> &iMaterialAssignment, const std::vector<cfgScalar> &iBaseMaterialDensities);
  cfgScalar computePoissonRatio(const std::vector<cfgScalar> &iStrainsAxialDir, const std::vector<cfgScalar> &iStrainsTransversalDir);
  cfgScalar computePoissonRatio(cfgScalar &iStrainsAxialDir, cfgScalar&iStrainsTransversalDir);
  cfgScalar computeVonMisesStress(std::vector<cfgScalar> &iStresses);  // iStresses = [sigma_xx sigma_yy sigma_zz sigma_xy sigma_yz sigma_zx]
  cfgScalar computeVonMisesStress2D(std::vector<cfgScalar> &iStresses); //iStresses = [sigma_xx sigma_yy sigma_xy]

  void fromCubicToOrthotropicProperties(int idim, const std::vector<cfgScalar> &iCubicProperties, std::vector<cfgScalar> &oOrthotropicProperties);

  // Hexahedral mesh utilities
  // -------------------------
  void getMaterialAssignment(int nx, int ny, int nz, const std::vector<int> &iMaterialCombIndices, int nX, int nY, int nZ, const std::vector<std::vector<int> > &iBaseCellMaterials, int repX, int repY, int repZ, int iNbSubdivisions, 
                             std::vector<int> &oMaterials);

  // Hexahedral mesh utilities 2D/3D
  // ----------------------------
  int getGridToVectorIndex(int i, int j, int nx, int ny);
  int getGridToVectorIndex(int i, int j, int k, int nx, int ny, int nz);
  void getVectorIndexToGrid(int iIndex, int nx, int ny, int &oIndex_i, int &oIndex_j);
  void getVectorIndexToGrid(int iIndex, int nx, int ny, int nz, int &oIndex_i, int &oIndex_j, int &oIndex_k);

  //0: Left, 1: Right, 2: Bottom, 3:Top
  void getSideElements(int iSide, int nx, int ny, std::vector<int> &oElementIndices);
  void getSideVertices(int iSide, int nx, int ny, std::vector<int> &oVertexIndices);
  void getSideVertices(int iSide, int nx, int ny, int nz, std::vector<int> &oVertexIndices);
  void getBoundaryCellIndices(int nx, int ny, std::vector<int> &oCellIndices);
  void getMaterialAssignment(int nx, int ny, const std::vector<int> &iMaterialCombIndices, int nX, int nY, const std::vector<std::vector<int> > &iBaseCellMaterials, int repX, int repY, int iNbSubdivisions, std::vector<int> &oMaterials);

  void repMaterialAssignment(int nx, int ny, const std::vector<int> &iMaterialCombIndices, int iNbSubdivions, int repX, int repY, std::vector<int> &oMaterials);
  void repMaterialAssignment(int nx, int ny, int nz, const std::vector<int> &iMaterialCombIndices, int repX, int repY, int repZ, int iNbSubdivisions, std::vector<int> &oMaterials);
  bool isVertexManifold(int indx, int indy, int nx, int ny,  const std::vector<int> &iMaterials);
  bool isVertexManifold(int indx, int indy, int indz, int nx, int ny,  int nz, const std::vector<int> &iMaterials, bool isMirrored, std::vector<int> *ioNonManifoldCells = NULL);
  bool isStructureManifold(int nx, int ny, const std::vector<int> &iMaterials, int nX, int nY, bool isStructuredMirrored, int iNbCellsToCheck);
  bool isStructureManifold(int nx, int ny, const std::vector<int> &iMaterials, bool isStructuredMirrored);
  bool isStructureManifold(int nx, int ny, int nz, const std::vector<int> &iMaterials, bool isStructuredMirrored, int nX, int nY, int nZ);
  bool isStructureManifold(int nx, int ny, int nz, const std::vector<int> &iMaterials, bool isStructuredMirrored, int nX, int nY, int nZ, int nmat);
  void getNonManifoldVertices(int nx, int ny, const std::vector<int> &iMaterials, bool isStructuredMirrored, std::vector<int> &oNonManifoldVertices);
  void getNonManifoldVertices(int nx, int ny, int nz, const std::vector<int> &iMaterials, bool isStructuredMirrored, std::vector<int> &oNonManifoldVertices);
  void updateMaterialSubstructure(std::vector<int> &ioMaterialAsssignment, int nx, int ny, const std::vector<int> &iSubMaterialStructure, int Nx, int Ny, int iSubMatStructLocation);
  
  void getLayer(int Nx, int NY, const std::vector<int> &iStructureElements, int iIndex, int iDim, std::vector<int> &oNewLayerElements);
  void getLayer(int Nx, int Ny, int Nz, const std::vector<int> &iStructureElements, int iIndex, int iDim, std::vector<int> &oNewLayerElements);
  void insertLayer(int Nx, int Ny, const std::vector<int> &iStructureElements, std::vector<int> &iNewLayerElements, int iIndex, int iDim, std::vector<int> &oNewStructureElements);
  void upscaleStructure(int Nx, int Ny, const std::vector<int> &iMatAssignment, std::vector<int> &oNewMaterialAssignments);
  void upscaleStructure(int Nx, int Ny, int Nz, const std::vector<int> &iMatAssignment, std::vector<int> &oNewMaterialAssignment);

  void mirrorStructure(int Nx, int Ny, const std::vector<int> &iStructureElements, std::vector<int> &oNewStructureElements);
  void mirrorStructure(int Nx, int Ny, int NZ, const std::vector<int> &iStructureElements, std::vector<int> &oNewStructureElements);
  void mirrorStructureAlongDiagonal(int Nx, int Ny, const std::vector<int> &iStructureElements, std::vector<int> &oNewStructureElements);
  void mirrorStructureAlongDiagonal(int Nx, int Ny, int Nz, const std::vector<int> &iStructureElements, std::vector<int> &oNewStructureElements);
  int getDiagonalStructureNumberOfElements(int N, int iDim);

  void getQuarter(int Nx, int Ny, const std::vector<int> &iStructureElements, std::vector<int> &oNewStructureElements);
  void getQuarter(int Nx, int Ny, int Nz, const std::vector<int> &iStructureElements, std::vector<int> &oNewStructureElements);
  void getTriangularStructure(int Nx, int Ny, const std::vector<int> &iStructureElements, std::vector<int> &oNewStructureElements);
  void getTetrahedralStructure(int Nx, int Ny, int Nz, const std::vector<int> &iStructureElements, std::vector<int> &oNewStructureElements);
  void getTriangleElements(int Nx, int Ny, const std::vector<int> &iStructureElements, std::vector<int> &oTriangleElements);
  void getTetrahedralElements(int Nx, int Ny, int Nz, const std::vector<int> &iStructureElements, std::vector<int> &oTetrahedralElements);

  void translateStructure(int Nx, int Ny, const std::vector<int> &iStructureElements, int Tx, int Ty, std::vector<int> &oTranslatedStructure);
  void translateStructure(int Nx, int Ny, int Nz, const std::vector<int> &iStructureElements, int Tx, int Ty, int Tz, std::vector<int> &oTranslatedStructure);

  cfgScalar computeStrain(const std::vector<cfgScalar> &ix, int nx, int ny, int iAxis);
  cfgScalar computeStrain(const std::vector<cfgScalar> &ix, int nx, int ny, int nz, int iAxis);
  void computeStrain(int n[2], const std::vector<std::vector<std::vector<cfgScalar> > > iX[2], std::vector<std::vector<cfgScalar> > oStrains[2][2]);
  void computeStrain3D(int n[3], const std::vector<std::vector<std::vector<cfgScalar> > > iX[3], std::vector<std::vector<cfgScalar> > oStrains[3][3]);
 
  int getMaterialSignature(int nx, int ny, const std::vector<int> &icellMaterials, const std::vector<int> &iIndicesToUse);
  void clusterByBoundary(int nx, int ny, const std::vector<std::vector<int> > &icellMaterials, int iNbMat, std::vector<std::vector<int> > &oClusters);

  void dumpStructure(int Nx, int Ny, const std::vector<int> &iMatAssignment);
  void dumpStructure(int Nx, int Ny, int Nz, const std::vector<int> &iMatAssignment);

  void getNeighbours(int Nx, int Ny, int indCell, std::vector<int> &oNeighbours, bool iUseTiling=true);
  void getDisconnectedComponents(int Nx, int Ny, const std::vector<int> &iMatAssignment, std::vector<std::vector<int> > &oComponents);
  bool isComponentConnectedToBorder(int Nx, int Ny, const std::vector<int> &iMatAssignment, std::vector<int> &iComponentcell);
  bool filterOutNonConnectedComponents(int Nx, int Ny, std::vector<int> &ioMatAssignment, bool &oModified);

  void getNeighbours(int Nx, int Ny, int Nz, int indCell, std::vector<int> &oNeighbours, bool iUseTiling=true);
  void getDisconnectedComponents(int Nx, int Ny, int Nz, const std::vector<int> &iMatAssignment, std::vector<std::vector<int> > &oComponents);
  bool isComponentConnectedToBorder(int Nx, int Ny, int Nz, const std::vector<int> &iMatAssignment, std::vector<int> &iComponentcell);
  bool filterOutNonConnectedComponents(int Nx, int Ny, int Nz, std::vector<int> &ioMatAssignment, bool &oModified);

  void convert2DCubicStructuresTo3DCubicStructures(int N, const std::vector<int> &iMatAssignment2D, std::vector<int> &oMatAssignment3D);

  // Point cloud utitlities
  // ----------------------
  void getFurthestPointsGreedy(int iOutputNbPoints, const std::vector<cfgScalar> &iPoints, int iDim, std::vector<int> &oPointIndices);

  void getBoundingBox(const std::vector<cfgScalar> &iPoints, int iDim, std::vector<cfgScalar> oBox[2]);
  bool rescaleData(std::vector<cfgScalar> &ioPoints, int iDim,  const std::vector<cfgScalar> &iTargetBoxLengths, std::vector<cfgScalar> *ioScalingFactors=NULL);
  void convertToLogValues(std::vector<cfgScalar> &ioPoints, int iDim, const std::vector<int> &iDimToScale, cfgScalar iEpsilon=1.e-6);

  // Triangle mesh utilities
  //------------------------
  float getMeanEdgeLength(const std::vector<float> &iX,  const std::vector<int> &iIndexArray, int iDim);
  Eigen::Vector3f getMeanEdgeLengthPerAxis(const std::vector<float> &iX, const std::vector<int> &iIndexArray, int iDim);
  Eigen::Vector3f getMedianEdgeLengthPerAxis(const std::vector<float> &iX, const std::vector<int> &iIndexArray, int iDim);
  float getMinEdgeLength(const std::vector<float> &iX,  const std::vector<int> &iIndexArray, int iDim);
  void getEdgesFromTriFaceIndexArray(const std::vector<int> &iTriIndexArray, std::vector<int> &oEdgeIndexArray);
  void sampleMesh(const std::vector<float> &iPoints, const std::vector<int> &iTriIndexArray, int iNbSamplesPerFace, std::vector<float> &oPoints);

  // Tet mesh utilities
  // ------------------
  void getEdgesFromTetFaceIndexArray(const std::vector<int> &iTetIndexArray, std::vector<int> &oEdgeIndexArray);

  // Conversion between material parameters
  void fromLamesParametersToYoungModulusPoissonRatio(cfgScalar iLambda, cfgScalar iMu, cfgScalar &oYoungModulus, cfgScalar &oPoissonRatio);
  void fromYoungModulusPoissonRatioToLamesParameters(cfgScalar iYoungModulus, cfgScalar iPoissonRatio, cfgScalar &oLambda, cfgScalar &oMu);

  // Conversion between types
  // -----------------------
  Eigen::Vector3f getVector3f(int indVertex, const std::vector<float> &iPoints);
  Eigen::Vector2f getVector2f(int indVertex, const std::vector<float> &iPoints);
  Vector3S getVector3S(int indVertex, const std::vector<cfgScalar> &iPoints);
  Vector2S getVector2S(int indVertex, const std::vector<cfgScalar> &iPoints);
  Vector3d getVector3d(int indVertex, const std::vector<double> &iPoints);
  Vector2d getVector2d(int indVertex, const std::vector<double> &iPoints);

  std::vector<float> toVectorFloat(const std::vector<Eigen::Vector3f> &iPoints);
  std::vector<float> toVectorFloat(const std::vector<Eigen::Vector2f> &iPoints);
  MatrixXS toMatrixXS(const std::vector<cfgScalar> &iValues);

  std::vector<cfgScalar> toVectorScalar(const std::vector<Vector3S> &iPoints);
  std::vector<cfgScalar> toVectorScalar(const std::vector<Vector2S> &iPoints);

  std::vector<Eigen::Vector2f> toVector2f(const std::vector<float> &iPoints);
  std::vector<Vector2S> toVector2S(const std::vector<cfgScalar> &iPoints);
  std::vector<Eigen::Vector3f> toVector3f(const std::vector<float> &iPoints);
   std::vector<Vector3S> toVector3S(const std::vector<cfgScalar> &iPoints);

  MatrixXS toMatrixScalar(const Eigen::MatrixXf &iMatrix);
  MatrixXS toMatrixScalar(const Eigen::MatrixXd &iMatrix);
  Eigen::MatrixXd toMatrixDouble(const MatrixXS &iMatrix);

  void toPardisoMatrix(Eigen::SparseMatrix<cfgScalar> &iMatrix, bool iTriangular, std::vector<int> &oRowIndicesRowMajor, std::vector<int> &oColIndicesRowMajor, std::vector<double> &oValues);

  // Misc
  // ----
  std::vector<int> genIncrementalSequence(int iMin, int iMax, int iStep=1);
  void sortValues(const std::vector<cfgScalar> &iValues, std::vector<int> &oOrderedIndices);
  void getEigenValues2x2SymmetricMatrix(const cfgScalar *iMatValues, cfgScalar oEigenValues[2]); //iMat entries should be in this order: xx, yy, xy
  void getEigenValues3x3SymmetricMatrix(const cfgScalar *iMatValues, cfgScalar oEigenValues[3]); //iMat entries should be in this order: xx, yy, zz, yz, xz, xy
};

#endif 

