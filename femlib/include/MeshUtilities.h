#ifndef MeshUtilities_h
#define MeshUtilities_h

#include <vector>
#include <map>
#include <Eigen/Dense>
#include "cfgDefs.h"

class ElementRegGrid;
class ElementRegGrid2D;

namespace meshUtil
{
  // Hexahedral mesh utilities
  // -------------------------
  int getElementIndex(int nx, int ny, int nz, int ii , int jj, int kk);

  //0: Left, 1: Right, 2: Bottom, 3:Top, 4: Back, 5: Front
  void getSideVoxels(int iSide, int nx, int ny, int nz, std::vector<int> &oSideVoxelIndices);
  void getSideVertices(int iSide, const ElementRegGrid * iElementGrid, std::vector<int> &oElemIndices, std::vector<std::vector<int> > &oFaceVertexIndices);
  void getSideVertices(int iSide, const ElementRegGrid2D * iElementGrid, std::vector<int> &oVertexIndices);
  void getExternalVertices(ElementRegGrid * iElementGrid, std::vector<int> &oElemIndices, std::vector<std::vector<int> > &oFaceVertexIndices, std::vector<Vector3S> &oNormals);
  void getVertexValences(const ElementRegGrid * iElementGrid, const std::vector<int> &iElemIndices, const std::vector<std::vector<int> > &iFvIndices, std::map<int,int> &oind2Valence);

  // Hexahedral mesh utilities 2D
  // ----------------------------
  //0: Left, 1: Right, 2: Bottom, 3:Top
  void getSideVertices(int iSide, const ElementRegGrid2D * iElementGrid, std::vector<int> &oElemIndices, std::vector<std::vector<int> > &oFaceVertexIndices);
  void getVertexValences(const ElementRegGrid2D * iElementGrid, const std::vector<int> &iElemIndices, const std::vector<std::vector<int> > &iEvIndices, std::map<int,int> &oind2Valence);
  cfgScalar getVonMisesStress(int iElementIndex, ElementRegGrid2D &iElementGrid);
  std::vector<cfgScalar> getVonMisesStressesPerElement(ElementRegGrid2D &iElementGrid);

  void computeCoarsenedElasticityTensor(ElementRegGrid2D &iElementGrid, const std::vector<std::vector<cfgScalar> > &iHarmonicDisplacements, MatrixXS &oCoarsenedTensor);
};

#endif 
