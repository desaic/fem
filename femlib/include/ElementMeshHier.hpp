#ifndef ELEMENTMESHHIER_HPP
#define ELEMENTMESHHIER_HPP
#include <vector>
#include "vecmath.h"
#include "BoundVec3.hpp"

class Element;
class Material;
struct MatrixXd;
class ElementMesh;

///@brief quantities stored in a quadrature point
struct QuadPt{
  ///@brief deformation gradient at each level
  std::vector<Matrix3f> F;
  ///@brief quadrature weight at each level
  std::vector<float> w;
  ///@brief natural coordinate at each level
  std::vector<Vector3f> X;
  ///@brief element indices
  std::vector<int> ei;
};

///@brief Assumes materials are MaterialQuad. 
///For each element, all quadrature points uses the first material stored in MaterialQuad.
class ElementMeshHier{
public:
  ///@brief Each mesh correspond to a level. 0 is the fine level.
  ///Materials, constrains and forces from the fine level are used.
  std::vector<ElementMesh*>m;
  ///@brief a list of quadrature points;
  std::vector<QuadPt> quadpt;
  ElementMeshHier();
  ///@brief collect quadrature points for one element
  std::vector<int> getQuadIdx(int level, int eIdx);
  ///@brief construct nLevel additional levels given m[0].
  void buildHier(int nLevel);

  ///@param qidx index of quadrature point
  void updateDefGrad(int level, int qidx);
  void updateDefGrad(int level);
  float getEnergy();
  float getEnergy(int level, int eIdx);

  ///@brief computes internal forces only
  std::vector<Vector3f> getForce(int level, int eIdx);
  std::vector<Vector3f> getForce(int level);

  virtual ~ElementMeshHier();
};

#endif
