#ifndef ELEMENTMESHHIER_HPP
#define ELEMENTMESHHIER_HPP
#include <vector>
#include "vecmath.h"
#include "BoundVec3.hpp"

class Element;
class Material;
struct MatrixXd;
class ElementMesh;

class ElementMeshHier{
public:
  ///@brief Each mesh correspond to a level. 0 is the fine level.
  ///Materials, constrains and forces from the fine level are used.
  std::vector<ElementMesh*>m;
  
  ElementMeshHier();

  ///@brief construct nLevel additional levels given m[0].
  void buildHier(int nLevel);

  float getEnergy();
  float getEnergy(int level, int eIdx);

  ///@brief computes internal forces only
  std::vector<Vector3f> getForce(int level, int eIdx);
  std::vector<Vector3f> getForce(int level);

  virtual ~ElementMeshHier();
};

#endif
